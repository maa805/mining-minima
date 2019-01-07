import numpy as np
import scipy
import itertools

from pyrosetta import *
from pyrosetta.rosetta import *
init()

class MiningMinima:
    '''
    A single MM style simulation
    '''

    def __init__(self, seq1='', seq2='', infile='', scorefxn = 'stepwise/rna/turner')

    '''
    Create a MM object.
    
    Parameters
    ---------------------------
    seq1:   The sequence for a single strand simulation or the first strand in a double strand simulation
    seq2:   The second sequence for a double strand simulation
    infile: A PDB file on which to run mining minima
    '''
        
		KT_IN_KCAL = 0.6163
		self.kT = 1.0
		
		if !infile: 

            self.seq1 = seq1
            if seq2: self.seq2 = seq2

        else: self.infile = infile
    
		self.input_pose = Pose()
		if seq1: self.movemap, self.dof_dict = pose_setup_turner(self, self.input_pose, seq1, seq2)
		else: self.movemap, self.dof_dict = pose_setup_from_file(self, self.input_pose, infile)
		
		N_DOFS = len(dof_dict)
		
		self.scorefxn = core.scoring.ScoreFunctionFactory.create_score_function(scorefxn)
		
		self.min_pose = Pose()
		find_minimum(self, self.min_pose, self.scorefxn)
		self.min_energy = self.scorefxn(self.min_pose)
		
		# Calculate hessian at base of minimum and normal modes
		self.hessian = calc_hessian_at_min(self, self.min_pose, self.scorefxn, self.dof_dict)
		self.eigenvalues, self.modes = np.linalg.eigh(self.hessian)
		
		# Calculate free energy
		self.harmonic_free_energy = self.min_energy - 0.5*N_DOFS*np.log(2.*np.pi*self.kT) + 0.5*np.log(np.product(self.eigenvalues))
		self.free_energy = mode_scan(self, self.min_pose, self.scorefxn, self.dof_dict)
		
	
	
	
	def pose_setup_turner(self, self.input_pose, self.seq1, self.seq2)
	
    def pose_setup_from_file(self, self.input_pose, self.infile) 
	
	def find_minimum(self, self.min_pose, self.scorefxn)
	
	def calc_hessian_at_min(self, self.min_pose, self.scorefxn, self.dof_dict)
	
	def mode_scan(self, self.min_pose, self.scorefxn, self.dof_dict)
	
	def calc_harmonic_free_energy
