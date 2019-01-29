from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.id import *
init()
# If necesarry, pose_setup_turner can revert to pose_setup_turner_ss and remove legacy to get parent 
# function

def add_bb_suite(suite, idx, movemap, dof_dict):

	dof_dict.update({
					5*idx: TorsionID(suite-1, BB, 5),5*idx+1: TorsionID(suite-1, BB, 6),
					5*idx+2: TorsionID(suite, BB, 1),5*idx+3: TorsionID(suite, BB, 2),  
					5*idx+4: TorsionID(suite, BB, 3)})
	
	for key in dof_dict: movemap.set(dof_dict[key], True)
	
	return movemap, dof_dict
	
def add_chi_dofs(n_residues, idx,  movemap, dof_dict):

	for ii in range(n_residues):
	
		dof_dict.update({5*(idx) + ii: TorsionID(ii + 1, CHI, 1)})
		movemap.set(TorsionID(ii + 1, CHI, 1), True)
		
	return movemap, dof_dict
	
def add_rigid_body_dofs(movemap, dof_dict):

	dof_dict.update({
					-1: numeric.xyz_Vector_double_t(1,0,0), -2: numeric.xyz_Vector_double_t(0,1,0},
					-3: numeric.xyz_Vector_double_t(0,0,1), -4: numeric.xyz_Vector_double_t(1,0,0),
					-5: numeric.xyz_Vector_double_t(0,1,0), -6: numeric.xyz_Vector_double_t(0,0,1)})
	
	movemap.set_jump(1, True)
	
	return 