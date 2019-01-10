from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.id import *
init()
# If necesarry, pose_setup_turner can revert to pose_setup_turner_ss and remove legacy to get parent 
# function
def pose_setup_turner_legacy(seq1, seq2):

    if seq2:
		if len(seq1) == len(seq2): pose, movemap, dof_dict = pose_setup_turner_ds(seq1, seq2)

		else: pose, movemap, dof_dict = pose_setup_turner_dangle(seq1, seq2)
	
	else: pose, movemap, dof_dict = pose_setup_turner_ss(seq1)
    return pose, movemap, dof_dict

def pose_setup_turner(seq):

	n_residues = len(seq1 + seq2)
	
	dof_dict = {}
	movemap = MoveMap()
	
	pose = protocols.recces.pose_setup_turner(seq1, seq2)
	
	res_idx = 0
	
	while res_idx < n_residues - 1:
	
		res_idx += 1
		suite = res_idx + 1
		
		if pose.fold_tree().is_cutpoint(res_idx):
			
			suite -= 1
			continue 
		
		movemap, dof_dict = add_bb_suite(suite, movemap, dof_dict)				
	
	movemap, dof_dict = add_chi_dofs(n_residues, movemap, dof_dict)
	
    return pose, movemap, dof_dict

def pose_setup_turner_ds(seq1, seq2):

	n_residues = len(seq1 + seq2)
	
	dof_dict = {}
	
	pose = protocols.recces.pose_setup_turner(seq1, seq2)
	
	for ii in range(n_residues - 1):
	
		if pose.fold_tree().is_cutpoint(ii+1):
		
		dof_dict.update({5*ii: TorsionID(ii+1, BB, 5)}, {5*ii+1: TorsionID(ii+1, BB, 6)
						{5*ii+2: TorsionID(ii+2, BB, 1). 
    return pose, movemap, dof_dict

def add_bb_suite(suite, movemap, dof_dict):
	
	dof_dict.update({5*suite: TorsionID(suite-1, BB, 5)}, {5*suite+1: TorsionID(suite-1, BB, 6),
					{5*suite+2: TorsionID(suite, BB, 1), {5*suite+3: TorsionID(suite, BB, 2)},  
					{5*suite+4: TorsionID(suite, BB, 3)})
	
	for key in dof_dict: movemap.set(dof_dict[key], True)
	
	return movemap, dof_dict
	
def add_chi_dofs(n_residues, movemap, dof_dict)

	for ii in range(n_residues):
	
		dof_dict.update({5*(n_residues - 1) + ii: TorsionID(ii + 1, CHI, 1)})
		movemap.set(TorsionID(ii + 1, CHI, 1))
		
	return movemap, dof_dict