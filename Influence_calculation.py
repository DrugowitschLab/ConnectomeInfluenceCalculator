import pandas as pd
import sqlite3
import numpy as np
from petsc4py import PETSc
from slepc4py import SLEPc

class InfluenceCalculator:
	def __init__(self, filename, signed=False):
		"""
		1. Load SQL data
		2. Create neuron_id -> W_id mapping, W_id -> neuron_id mapping (using bidict)
		3. Create sparse W (unnormalized, eventually signed)
		"""
		self.signed_W = signed
		pass

	def _load_sql_data(self, filename):

		# Connect to the SQLite database
		conn = sqlite3.connect(filename)

		# List all tables in the database
		cursor = conn.cursor()
		cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
		tables = cursor.fetchall()

		# Print the tables
		print("Tables in the database:", tables)

		# Get the meta data, cell types, etc.
		meta = pd.read_sql_query("SELECT * FROM meta", conn)

		### Construct the SQL query for the edgelist and add condition on minimum synaptic count (here, min=5)
		query = f"""
		SELECT *
		FROM edgelist_simple
		WHERE count >= 5
		"""

		# Execute the query and collect the results
		elist = pd.read_sql_query(query, conn)

		# Add 'post_count' column to elist
		elist['post_count'] = np.round(elist['count']/elist['norm'])

		# Close the database connection
		conn.close()
		
		return meta, elist

	def _create_neuron_W_id_mapping(self, elist):

		# Find unique neuron ids
		unique_ids = pd.unique(np.hstack([elist['post'].to_numpy(), elist['pre'].to_numpy()]))
		# Number of neurons
		n_neurons = len(unique_ids)

		# Create a mapping between matrix indices and neuron ids
		id_to_index = {id_: index for index, id_ in enumerate(unique_ids)}
		index_to_id = {index: id_ for id_, index in id_to_index.items()}

		self.id_to_index = id_to_index
		self.index_to_id = index_to_id
		self.n_neurons = n_neurons

	def _create_sparse_W(self, meta, elist, syn_weight_measure='norm'):

		# If signed W, change synaptic count to negative in elist for inhibitory neurons
		if self.signed_W:
			negative_nt = {'glutamate', 'gaba', 'serotonin', 'octopamine'}
			# Create a boolean mask for rows in meta that meet our conditions
			mask = meta['top_nt'].isin(negative_nt) & meta['id'].isin(elist['pre'])
			# Get the ids that need to be updated
			ids_to_update = set(meta.loc[mask, 'id'])
			# Update elist in one vectorized operation
			elist.loc[elist['pre'].isin(ids_to_update), 'count'] *= -1

		# Get synaptic weights
		syn_weights = elist[syn_weight_measure].values

		# Map pre and post indices
		pre_ind = elist['pre'].map(self.id_to_index)
		post_ind = elist['post'].map(self.id_to_index)

		# Remove connections involving neurons not in the specified dataset
		mask = (~pre_ind.isna()) & (~post_ind.isna())
		pre_ind = pre_ind[mask].astype(int)
		post_ind = post_ind[mask].astype(int)
		syn_weights = syn_weights[mask]

		"""Create a sparse random matrix using PETSc"""
		matrix = PETSc.Mat().create()
		matrix.setSizes([self.n_neurons, self.n_neurons])
		matrix.setFromOptions()
		matrix.setType('aij')  # sparse matrix
		matrix.setUp()

		for i, j, v in zip(post_ind, pre_ind, syn_weights):
			matrix.setValue(i, j, v, addv=True)

		matrix.assemblyBegin()
		matrix.assemblyEnd()
		
		self.W = matrix

	def _normalize_W(W_norm):

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

		# Create W = alpha*connectivity_matrix - I
		alpha = 0.99 / eig_val_largest
		W_norm.scale(alpha)
		W_norm.shift(-1.0)

		return W_norm

	def _solve_lin_system(W_norm, s):
		"""
		Solves the system W_norm * x = s
		"""
		
		b = PETSc.Vec().createWithArray(s, comm=PETSc.COMM_SELF)

		ksp = PETSc.KSP().create()
		ksp.setOperators(W_norm)
		ksp.setType(PETSc.KSP.Type.GMRES)  # or another appropriate solver type
		pc = ksp.getPC()
		pc.setType(PETSc.PC.Type.ILU)  # or another appropriate preconditioner

		x = W_norm.createVecRight()
		ksp.solve(b, x)

		return x
	
	def _set_columns_to_zero(self, silenced_neurons):
		# Copy W
		W_norm = self.W.copy()

		# Get matrix dimensions
		m, n = W_norm.getSize()
		
		# Convert silenced_neurons to 32-bit integers and ensure it's a NumPy array
		silenced_neurons_32 = np.asarray(silenced_neurons, dtype=np.int32)
		
		# Create a scaling vector initialized with ones
		scale_vec = W_norm.createVecRight()
		scale_vec.set(1.0)
		
		# Set the entries corresponding to silenced_neurons to zero in the scaling vector
		scale_vec.setValues(silenced_neurons_32, np.zeros_like(silenced_neurons_32))
		scale_vec.assemble()
		
		# Scale the matrix columns
		W_norm.diagonalScale(None, scale_vec)
		
		# Assemble the matrix after modifications
		W_norm.assemble()
		
		return W_norm

	def calculate_influence(self, seed_ids, silenced_neurons=[]):

		# map seed_ids to W_ids to get seed_vec		
		seed_vec = np.zeros(self.n_neurons)
		seed_indices = []
		for ii in seed_ids:
			if ii in self.id_to_index:
				seed_indices.append(self.id_to_index[ii])
		seed_vec[seed_indices] = 1

		# Inhibit specific neurons
		if silenced_neurons:
			# Map to W indices and exclude seed neurons 
			silenced_indices_temp = np.array([self.id_to_index[id] for id in silenced_neurons if id in self.id_to_index])
			exclusion_indices = np.where(seed_vec==1)[0]
			silenced_indices = np.setdiff1d(silenced_indices_temp, exclusion_indices)
			W_norm = self._set_columns_to_zero(self, silenced_indices)

		# normalize W_norm and remove I
		W_norm = self._normalize_W(W_norm)

		influence_vec = self._solve_lin_system(W_norm, -seed_vec)

		# turn influence_vec into panda table with rows (neuron_id, influence, part_of_seed)



ic = InfluenceCalculator('filename')
# load sql file
# create neuro_id -> W id mapping
# create connectivity matrix W + normalize