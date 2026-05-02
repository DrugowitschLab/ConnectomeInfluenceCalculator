import pandas as pd
import sqlite3
import numpy as np
from petsc4py import PETSc
from slepc4py import SLEPc
from bidict import bidict


class InfluenceCalculator:
    def __init__(self, filename, signed=False, count_thresh=5,
                 syn_weight_measure='count', inhibitory_nts=None,
                 excluded_nts=None, lambda_max=0.99):
        """
        Creates a class instance by loading SQL data with the given
        filename, establishing a neuron_id <-> W_id mapping (using
        bidict), and by creating a sparse W matrix that contains the
        synaptic weights drawn from the column named by
        syn_weight_measure (signed if signed=True).

        syn_weight_measure selects which edge column populates W: 'count'
        is the raw synapse count, 'norm' is the per-postsynaptic input
        fraction (count / sum(count) per post). The default is 'count' so
        that the signed=True negation has a clean interpretation; see
        signed below. No silent default of 'norm' is provided -- callers
        are expected to choose deliberately.

        signed=True multiplies the chosen syn_weight_measure column by -1
        for edges whose pre-neuron's top_nt is in inhibitory_nts. Note
        that flipping the sign of 'norm' values means the columns of W no
        longer sum to 1, so the input-normalisation interpretation is
        lost; 'count' is the more natural choice when signed=True.

        inhibitory_nts is a set or list of neurotransmitter names
        (matching values in the 'top_nt' metadata column) that should be
        treated as inhibitory. Required when signed=True; ignored when
        signed=False. The library has no per-organism default -- the
        caller must supply the set explicitly (e.g. {'gaba'} for
        C. elegans, or {'glutamate', 'gaba', 'serotonin', 'octopamine'}
        for the historical Drosophila convention).

        excluded_nts is a set or list of neurotransmitter names whose
        pre-neurons contribute nothing to W: their outgoing edges are
        removed from the connectivity matrix entirely. Independent of
        signed=True/False. Use this for transmitter classes whose net
        sign at a given target depends on the receptor mix and so cannot
        be assigned a single sign safely (e.g. dopamine, serotonin,
        octopamine in C. elegans).

        lambda_max is the target largest real eigenvalue of the rescaled
        W after normalisation; W is scaled in place by
        lambda_max / lambda_max(W) so that lambda_max of the rescaled W
        equals lambda_max exactly. Must satisfy 0 < lambda_max < 1 for
        the steady-state solve to remain stable. The amplification of the
        leading eigenmode in (I - W_rescaled)^-1 is 1 / (1 - lambda_max),
        so the default 0.99 gives ~100x and a smaller value
        (e.g. 0.5 -> 2x) damps the global mode and exposes per-target
        seed-specificity at the cost of attenuating long polysynaptic
        paths.
        """
        if signed and inhibitory_nts is None:
            raise ValueError(
                "signed=True requires inhibitory_nts to be specified "
                "as a set of neurotransmitter names matching values in "
                "the 'top_nt' metadata column."
            )
        if not (0 < lambda_max < 1):
            raise ValueError(
                f"lambda_max must satisfy 0 < lambda_max < 1; got "
                f"{lambda_max}."
            )
        if syn_weight_measure not in ('count', 'norm'):
            raise ValueError(
                f"syn_weight_measure must be 'count' or 'norm'; got "
                f"{syn_weight_measure!r}."
            )

        self.W_signed = signed
        self.lambda_max = lambda_max
        self.syn_weight_measure = syn_weight_measure
        self.inhibitory_nts = (set(inhibitory_nts)
                               if inhibitory_nts is not None else set())
        self.excluded_nts = (set(excluded_nts)
                             if excluded_nts is not None else set())

        elist = self._load_sql_data(filename, count_thresh)
        self._create_neuron_W_id_mapping(elist)
        self._create_sparse_W(elist)

    def _load_sql_data(self, filename, count_thresh):
        """This method opens an SQLite database with the given filename,
        loads and stores metadata, and loads and returns a list of edges
        in the connectivity graph.
        """
        # Connect to the SQLite database
        conn = sqlite3.connect(filename)

        # List all tables in the database
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()

        # Print the tables
        # print("Tables in the database:", tables)

        # Get the meta data, cell types, etc.
        self.meta = pd.read_sql_query("SELECT * FROM meta", conn)

        # Construct the SQL query for the edgelist and add condition on
        # minimum synaptic count (here, min=count_thresh)
        query = f"""
        SELECT *
        FROM edgelist_simple
        WHERE count >= {count_thresh}
        """

        # Execute the query, collect the results, and close the db connection
        elist = pd.read_sql_query(query, conn)
        conn.close()

        # Add 'post_count' column to elist and return
        elist['post_count'] = np.round(elist['count']/elist['norm'])
        return elist

    def _create_neuron_W_id_mapping(self, elist):
        """This method uses the list of edges to find unique neuron ids,
        and create a bidirectional mapping from neuron IDs to rows and
        columns in the W matrix.
        """
        # Find unique neuron ids
        unique_ids = pd.unique(np.hstack([elist['post'].to_numpy(),
                                          elist['pre'].to_numpy()]))
        # Number of neurons
        self.n_neurons = len(unique_ids)

        # Create a bidirectional mapping from neuron IDs to matrix indices.
        self.id_to_index = bidict(
            {neuron_id: idx for idx, neuron_id in enumerate(unique_ids)})

    def _create_sparse_W(self, elist):
        """This method takes the edge list and uses it to populate the
        sparse connectivity matrix W from the column named by
        self.syn_weight_measure ('count' or 'norm').
        """
        syn_weight_measure = self.syn_weight_measure

        # Drop edges originating from neurons whose top_nt is in
        # excluded_nts; these contribute nothing to W regardless of
        # signed=True/False.
        if self.excluded_nts:
            if 'top_nt' not in self.meta.columns:
                raise ValueError(
                    "excluded_nts requires the SQLite 'meta' table to "
                    "include a 'top_nt' column identifying "
                    "neurotransmitter types."
                )
            mask = self.meta['top_nt'].isin(self.excluded_nts)
            excl_ids = set(self.meta.loc[mask, 'root_id'])
            elist = elist[~elist['pre'].isin(excl_ids)].copy()

        # If W ought to be signed, negate the weights for edges from
        # inhibitory pre-neurons.  We negate the syn_weight_measure
        # column actually consumed below; the original implementation
        # always negated 'count' regardless of which column populated W,
        # so the signed flag had no effect when 'norm' was used to build
        # W.  Note that when syn_weight_measure='norm' the negated
        # entries no longer leave the postsynaptic input fractions
        # summing to 1, so the input-normalisation interpretation is
        # lost; 'count' is the more natural choice when signed=True.
        if self.W_signed:
            if 'top_nt' not in self.meta.columns:
                raise ValueError(
                    "signed=True requires the SQLite 'meta' table to "
                    "include a 'top_nt' column identifying "
                    "neurotransmitter types."
                )
            mask = (self.meta['top_nt'].isin(self.inhibitory_nts) &
                    self.meta['root_id'].isin(elist['pre']))
            ids_to_update = set(self.meta.loc[mask, 'root_id'])
            elist.loc[elist['pre'].isin(ids_to_update),
                      syn_weight_measure] *= -1

        # Get synaptic weights
        syn_weights = elist[syn_weight_measure].values

        # Map pre and post indices
        pre_ind = elist['pre'].map(self.id_to_index)
        post_ind = elist['post'].map(self.id_to_index)

        # Create a sparse matrix using PETSc, and populate from edge list
        W = PETSc.Mat().create()
        W.setSizes([self.n_neurons, self.n_neurons])
        W.setFromOptions()
        W.setType('aij')  # sparse matrix
        W.setUp()

        for i, j, v in zip(post_ind, pre_ind, syn_weights):
            W.setValue(i, j, v, addv=True)

        W.assemblyBegin()
        W.assemblyEnd()

        self.W = W

    def _normalize_W(self, W_norm):
        """Rescale W_norm matrix so its largest real eigenvalue equals
        self.lambda_max, then subtract the identity matrix.  This ensures
        that the largest real eigenvalue of W_norm becomes
        lambda_max - 1 (negative for any 0 < lambda_max < 1) so the
        steady-state solve is stable.

        Always rescales (rather than only capping when the natural
        eigenvalue exceeds the target), so lambda_max is a true control
        knob over leading-mode amplification rather than just a stability
        ceiling.

        Please note that this method directly modifies W_norm without
        creating a copy.
        """
        # Compute the largest eigenvalue
        eps = SLEPc.EPS().create()
        eps.setOperators(W_norm)
        eps.setProblemType(SLEPc.EPS.ProblemType.NHEP)
        eps.setDimensions(1)
        eps.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
        eps.setFromOptions()
        eps.solve()

        if eps.getConverged() < 1:
            raise RuntimeError("Eigenvalue solver did not converge")

        eig_val_largest = eps.getEigenvalue(0).real

        # Rescale so the largest real eigenvalue equals self.lambda_max
        # exactly, then form W = alpha * W - I.
        if eig_val_largest > 0:
            alpha = self.lambda_max / eig_val_largest
            W_norm.scale(alpha)
        W_norm.shift(-1.0)

    @staticmethod
    def _solve_lin_system(W_norm, s):
        """Solves the system W_norm * x = s and returns x.
        """
        b = PETSc.Vec().createWithArray(s, comm=PETSc.COMM_SELF)

        ksp = PETSc.KSP().create()
        ksp.setOperators(W_norm)
        ksp.setType(PETSc.KSP.Type.GMRES)  # could use another solver type
        pc = ksp.getPC()
        pc.setType(PETSc.PC.Type.ILU)  # could use another preconditioner

        x = W_norm.createVecRight()
        ksp.solve(b, x)

        return x

    def _set_columns_to_zero(self, silenced_W_idcs):
        """This method returns a copy of W with the columns listed in
        the silenced_W_idcs (given by W indices of silenced neurons) set
        to zero.
        """
        # Copy W
        W_norm = self.W.copy()

        # Get matrix dimensions
        m, n = W_norm.getSize()

        # Convert silenced_W_idcs to 32-bit ints and ensure it's a NumPy array
        silenced_W_idcs_32 = np.asarray(silenced_W_idcs, dtype=np.int32)

        # Create a scaling vector initialized with ones
        scale_vec = W_norm.createVecRight()
        scale_vec.set(1.0)

        # Set silenced_W_idcs entries to zero in the scaling vector
        scale_vec.setValues(silenced_W_idcs_32,
                            np.zeros_like(silenced_W_idcs_32))
        scale_vec.assemble()

        # Scale the matrix columns
        W_norm.diagonalScale(None, scale_vec)

        # Assemble the matrix after modifications
        W_norm.assemble()

        return W_norm

    def _build_influence_dataframe(self, influence_vec, seed_vec):
        """This method turns the influence vector influence_vec into a
        pandas dataframe that lists neuron_ids with influence score and
        whether each neuron has been a seed neuron.  In signed mode the
        real part is preserved so that net-inhibited targets carry a
        negative score; in unsigned mode the magnitude is taken.
        """
        if self.W_signed:
            influence_vec = np.real(influence_vec)
        else:
            influence_vec = np.abs(np.real(influence_vec))

        seed_indices = np.where(seed_vec == 1)[0]
        seed_ids = np.array([self.id_to_index.inv[id] for id in seed_indices
                             if id in self.id_to_index.inv])

        # Build dataframe
        influence_df = pd.DataFrame({
            'matrix_index': list(self.id_to_index.inv.keys()),
            'id': list(self.id_to_index.inv.values()),
        })
        # Add is_seed column
        influence_df['is_seed'] = False
        influence_df.loc[influence_df['id'].isin(seed_ids), 'is_seed'] = True

        # Add the new influence score column
        W_signed_str = 'signed' if self.W_signed else 'unsigned'
        column_name = f"Influence_score_({W_signed_str})"
        influence_df[column_name] = influence_vec

        return influence_df

    def calculate_influence(self, seed_ids, silenced_neurons=[]):
        """This method calculates the influence score for the given
        seed id and list of silenced neurons. It returns the results
        as a pandas dataframe.
        """
        # map seed_ids to W_ids to get seed_vec
        seed_vec = np.zeros(self.n_neurons)
        seed_indices = []
        for ii in seed_ids:
            if ii in self.id_to_index:
                seed_indices.append(self.id_to_index[ii])
        seed_vec[seed_indices] = 1

        # Inhibit specific neurons
        if len(silenced_neurons) > 0:
            # Map to W indices and exclude seed neurons
            silenced_indices_temp = np.array(
                [self.id_to_index[id] for id in silenced_neurons
                 if id in self.id_to_index])
            exclusion_indices = np.where(seed_vec == 1)[0]
            silenced_indices = np.setdiff1d(silenced_indices_temp,
                                            exclusion_indices)
            W_norm = self._set_columns_to_zero(silenced_indices)
        else:
            # If no silencing
            W_norm = self.W.copy()

        self._normalize_W(W_norm)
        influence_vec = self._solve_lin_system(W_norm, -seed_vec)
        influence_df = self._build_influence_dataframe(influence_vec, seed_vec)

        return influence_df
