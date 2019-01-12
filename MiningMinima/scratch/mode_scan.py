import numpy as np
import scipy
from pyrosetta import *
from pyrosetta.rosetta import *
init()
def mode_scan(min_pose, mode, scorefxn, dof_dict, h = 1.0, range = 60.0, kt = 1.0):
	'''
	Given an input pose (assumed to be a local minimum) and a correspondig mode, returns the partition
	function evaluated along the mode.'''

## Initialize dofs from input pose and dof_dict -- list comprehension?
## We should probably return the scanned modes for manipulation, as well as a calculated free energy

	# Initialize stuff 
	min_energy = scorefxn(min_pose)
	energy = min_energy
	
	temp = Pose()
	temp.assign(min_pose)
	
	# Read off dofs from input pose
	dofs = [min_pose.torsion(dof_dict[key]) for key in dof_dict]
	
	ii = 0
	result = np.array([])
	
	# h is in degrees to match rosetta movers
	h = np.arange(-range, range+h, h)
	
	# Scan along mode at indicated granularity 
	for hh in h:
	
		dofs_new = dofs + hh*mode
		for key, val in enumerate(dofs_new): temp.set_torsion(dof_dict[key], val)
		
		result = np.append(result, np.exp(-(scorefxn(temp) - min_energy))/kt)
	
	# Apply trapezoidal rule to get anharmonic free energy
	partition_function = np.trapz(result, x=h*np.pi/180.0)
	
	return partition_function, result 
	
	
	