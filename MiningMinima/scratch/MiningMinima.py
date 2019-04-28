import numpy as np
import scipy
import itertools

from pyrosetta import *
from pyrosetta.rosetta import *
init()

from pose_setup import add_bb_suite
from pose_setup import add_chi_dofs
from calc_hessian_at_min import *
from mode_scan import *

class MiningMinima:	
	'''
	A single MM style simulation
	'''

	def __init__(self, seq1='', seq2='', infile='', scorefxn = 'stepwise/rna/turner', kt = 1.0):
		'''
		Create a MM object.
	
		Parameters
		---------------------------
		seq1:   The sequence for a single strand simulation or the first strand in a double strand simulation
		seq2:   The second sequence for a double strand simulation
		infile: A PDB file on which to run mining minima
		'''
		self.kt = kt
		
		self.harmonic_free_energy = 0.0
		self.anharmonic_free_energy = 0.0
		KT_IN_KCAL = 0.6163
		
		self.input_pose = Pose()
		if not infile: 
			
			self.seq1 = seq1
			self.seq2 = seq2
			self.pose_setup_turner()
			
		else: 
			self.infile = infile
			self.pose_setup_from_file()
		
		self.n_dofs = len(self.dof_dict)
		
		self.scorefxn = core.scoring.ScoreFunctionFactory.create_score_function(scorefxn)
		
		self.min_pose = Pose()
		self.min_pose.assign(self.input_pose)
		
		self.find_minimum()
		self.min_energy = self.scorefxn(self.min_pose)
		
		# Calculate hessian at base of minimum and normal modes
		self.hessian = calc_hessian_at_min(self.min_pose, self.scorefxn, self.dof_dict)
		self.eigenvalues, self.modes = np.linalg.eigh(self.hessian)
		
		# Function for density of states
		self.dos = lambda E: (2.0*np.pi)**(self.n_dofs/2)*(E-self.min_energy)**(self.n_dofs/2 -1 )/scipy.special.gamma(self.n_dofs/2)/np.sqrt(np.linalg.det(self.hessian))*np.heaviside(E-self.min_energy, 0.5)
		#self.dos = lambda E: 0.5/self.n_dofs*(2.0*np.pi*(E-self.min_energy))**(self.n_dofs/2)/scipy.special.gamma(self.n_dofs/2)/np.sqrt(np.linalg.det(self.hessian))*np.heaviside(E-self.min_energy, 0.5)
		
		# Calculate free energy
		self.calc_harmonic_free_energy()
		self.calc_anharmonic_free_energy()
		
	
	

	def pose_setup_turner(self):
		
		n_residues = len(self.seq1)
		if self.seq2: n_residues += len(self.seq2)
		
		dof_dict = {}
		movemap = MoveMap()
		
		pose = protocols.recces.pose_setup_turner(self.seq1, self.seq2)
		
		idx = 0
		suite = 2
		
		for i in range(n_residues - 1):
			
			if pose.fold_tree().is_cutpoint(i+1):
				
				suite += 1
				
				continue 
			
			movemap, dof_dict = add_bb_suite(suite, idx, movemap, dof_dict)				
			idx += 1
			suite += 1
		
		print dof_dict
			
		
		movemap, dof_dict = add_chi_dofs(n_residues, idx, movemap, dof_dict)

		self.input_pose.assign(pose)
		self.movemap = movemap
		self.dof_dict = dof_dict
		
		
	
	def pose_setup_from_file(self): 
		
		pose = pose_from_file(self.infile)
		
		n_residues = len(pose.sequence())
		ft = pose.fold_tree()
		dof_dict = {}
		movemap = MoveMap()
		
		idx = 0
		suite = 2
		
		for i in range(n_residues - 1):
		
			if pose.fold_tree().is_cutpoint(i+1):
			
				suite += 1
				
				continue
			
			movemap, dof_dict = add_bb_suite(suite, idx, movemap, dof_dict)
			idx += 1
			suite += 1
			
		movemap, dof_dict = add_chi_dofs(n_residues, idx, movemap, dof_dict)
		
		self.input_pose.assign(pose)
		self.movemap = movemap
		self.dof_dict = dof_dict
		
		return
		
	#def calc_hessian_at_min(self, self.min_pose, self.scorefxn, self.dof_dict)
	
	
	#def mode_scan(self, self.min_pose, self.scorefxn, self.dof_dict) 
	
	def find_minimum(self):
		'''Minimizes an input pose using indicated scorefunction and movemap'''
		
		minmover = rosetta.protocols.minimization_packing.MinMover(
				self.movemap, self.scorefxn,'lbfgs_armijo_nonmonotone', 1.0e-7, True)
		minmover.max_iter(1000000)
		minmover.min_options().use_nblist(True)
		minmover.min_options().nblist_auto_update(True)
		
		minmover.apply(self.min_pose)
		
	def calc_anharmonic_free_energy(self):
		
		total_partition = 1.0
		
		for i in range(self.n_dofs): 
			
			mode_partition, _ = mode_scan(self.min_pose, self.modes[:,i], self.scorefxn, self.dof_dict)
			total_partition *= mode_partition
			
		self.anharmonic_free_energy = self.min_energy - self.kt*np.log(total_partition)

	def calc_harmonic_free_energy(self):
		
		self.harmonic_free_energy = (self.min_energy - 0.5*self.kt*self.n_dofs*np.log(2*np.pi*self.kt) + 0.5*self.kt*np.log(np.prod(self.eigenvalues)))
	
	def harmonic_ensemble(self, n_struct=200):
		'''
		
		Parameters
		--------------------
		n_struct:	Number of ensemble members to generate
		'''
		
		ensemble = np.zeros((self.n_dofs, n_struct))
		
		mu = [self.min_pose.torsion(self.dof_dict[key]) for key in self.dof_dict]
		cov = np.matmul( np.matmul( self.modes, np.diag(1/self.eigenvalues)), self.modes.T)*(180./np.pi)**2	
		
		ensemble = np.random.multivariate_normal(mu, cov, size=(n_struct))
		
		return ensemble
#def pose_setup_from_file(pose, infile)
#
#	#TODO: Figure out how to determine whether is ss, ds, or something else and process accordingly 
