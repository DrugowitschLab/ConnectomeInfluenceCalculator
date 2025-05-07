import pandas as pd
import sqlite3
import numpy as np
from petsc4py import PETSc
from slepc4py import SLEPc
from bidict import bidict


# list of neurotransmitters that will receive negative signs if requested
NEG_NEUROTRANSMITTERS = {'glutamate', 'gaba', 'serotonin', 'octopamine'}


class InfluenceCalculator:
    def __init__(self, filename, signed=False, count_thresh=5):
        """
        Creates a class instance by loading SQL data with the given
        filename, establishing a neuron_id <-> W_id mapping (using bidict), and
        by creating a sparse W matrix that contains the synapse count
        (signed if signed=True).
        """
        self.W_signed = signed
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

        # Construct the SQL query for the edgelist and add condition on minimum
        # synaptic count (here, min=5)
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

    def _create_sparse_W(self, elist, syn_weight_measure='norm'):
        """This method takes the edge list, and uses it to populate the
        sparse connectivity matrix W.
        syn_weight_measure takes either 'norm' for normalized postsynaptic 
        weights or 'count' for unnormalized postsynaptic weights
        """
        # If W ought to be signed, change relevant edge list entries
        if self.W_signed:
            # Create a boolean mask for rows in meta that meet our conditions
            mask = (self.meta['top_nt'].isin(NEG_NEUROTRANSMITTERS) &
                    self.meta['id'].isin(elist['pre']))
            # Get the ids that need to be updated
            ids_to_update = set(self.meta.loc[mask, 'id'])
            # Update elist in one vectorized operation
            elist.loc[elist['pre'].isin(ids_to_update), 'count'] *= -1

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
        """Rescale W_norm matrix to ensure that largest real eigenvalue is
        0.99, and then subtract identity matrix. This ensures that the largest
        real eigenvalue of W_norm is -0.01.

        Please note that this method directly modifies W_norm without creating
        a copy.
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

        # Create W = alpha * W - I
        alpha = 0.99 / eig_val_largest
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
        the silenced_W_idcs (given by W indices of silenced neurons) set to zero.
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
        """This method turns the influence vector influence_vec and turns it
        into a pandas dataframe that lists neuron_ids with influence score
        and whether each neuron has been a seed neuron.
        """
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
