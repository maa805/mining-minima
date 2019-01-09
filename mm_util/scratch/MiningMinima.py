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
        
		self.harmonic_free_energy = 0.0
		self.anharmonic_free_energy = 0.0
		KT_IN_KCAL = 0.6163
		
		if !infile: 

            self.seq1 = seq1
            if seq2: self.seq2 = seq2

        else: self.infile = infile
    
		self.input_pose = Pose()
		if seq1: self.movemap, self.dof_dict = pose_setup_turner(self, self.input_pose, seq1, seq2)
		else: self.movemap, self.dof_dict = pose_setup_from_file(self, self.input_pose, infile)
		
		self.n_dofs = len(dof_dict)
		
		self.scorefxn = core.scoring.ScoreFunctionFactory.create_score_function(scorefxn)
		
		self.min_pose = Pose()
		find_minimum(self, self.min_pose, self.scorefxn)
		self.min_energy = self.scorefxn(self.min_pose)
		
		# Calculate hessian at base of minimum and normal modes
		self.hessian = calc_hessian_at_min(self, self.min_pose, self.scorefxn, self.dof_dict)
		self.eigenvalues, self.modes = np.linalg.eigh(self.hessian)
		
		# Calculate free energy
		self.calc_harmonic_free_energy()
		self.calc_anharmoic_free_energy()
		
	
	
	
	def pose_setup_turner(self, self.input_pose, self.seq1, self.seq2)
	
		return movemap, dof_dict
	
    def pose_setup_from_file(self, self.input_pose, self.infile) 
	
		return movemap, dof_dict
		
	def find_minimum(self, self.min_pose, self.scorefxn)
	#def calc_hessian_at_min(self, self.min_pose, self.scorefxn, self.dof_dict)
	
	
	#def mode_scan(self, self.min_pose, self.scorefxn, self.dof_dict) 
	def calc_anharmonic_free_energy(self):
	
		total_partition = 1.0
		
		for i in range(self.n_dofs): 
			
			mode_partition, _ = mode_scan(self.min_pose, self.scorefxn, self.dof_dict)
			total_partition *= mode_partition
		
		self.anharmonic_free_energy = total_partition
			
		
	def calc_harmonic_free_energy(self):
	
		self.harmonic_free_energy = -kt*(self.min_energy - 0.5*self.n_dofs*np.log(2*np.pi*kt) + 0.5*self.n_dofs*np.log(np.prod(self.eigenvalues)))

def pose_setup_from_file(pose, infile)

	#TODO: Figure out how to determine whether is ss, ds, or something else and process accordingly 

def find_minimum(min_pose, scorefxn, movemap)

	#HMMMM maybe let user specify a particular minmover, otherwise just setup a default minmover? we'll worry about this later
		minmover = MinMover(
	
