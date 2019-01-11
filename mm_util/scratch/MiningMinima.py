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
		
		if not infile: 
			
			self.seq1 = seq1
			if seq2: self.seq2 = seq2
			
		else: self.infile = infile
		
		self.input_pose = Pose()
		self.pose_setup_turner()
		#else: self.movemap, self.dof_dict = pose_setup_from_file(self)
		
		self.n_dofs = len(self.dof_dict)
		
		self.scorefxn = core.scoring.ScoreFunctionFactory.create_score_function(scorefxn)
		
		self.min_pose = Pose()
		self.min_pose.assign(self.input_pose)
		
		self.find_minimum()
		self.min_energy = self.scorefxn(self.min_pose)
		
		# Calculate hessian at base of minimum and normal modes
		self.hessian = calc_hessian_at_min(self.min_pose, self.scorefxn, self.dof_dict)
		self.eigenvalues, self.modes = np.linalg.eigh(self.hessian)
		
		# Calculate free energy
		self.calc_harmonic_free_energy()
		self.calc_anharmonic_free_energy()
		
	
	

	def pose_setup_turner(self):
		
		n_residues = len(self.seq1 + self.seq2)
		
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
		print dof_dict
		self.input_pose.assign(pose)
		self.movemap = movemap
		self.dof_dict = dof_dict
		
		
	
	def pose_setup_from_file(self): 
		
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
		kt = 1.0
		self.harmonic_free_energy = (self.min_energy - 0.5*self.kt*self.n_dofs*np.log(2*np.pi*self.kt) + 0.5*self.kt*np.log(np.prod(self.eigenvalues)))

#def pose_setup_from_file(pose, infile)
#
#	#TODO: Figure out how to determine whether is ss, ds, or something else and process accordingly 
