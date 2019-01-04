from pyrosetta import *
from pyrosetta.rosetta import *
init()

def pose_setup_turner(seq1, seq2):

    if len(seq1) == len(seq2): pose, movemap, dof_dict = pose_setup_turner_ds(seq1, seq2)

    else: pose, movemap, dof_dict = pose_setup_turner_ss(seq1)

    return pose, movemap, dof_dict



def pose_setup_turner_ss(seq):

    return pose, movemap, dof_dict




def pose_setup_turner_ds(seq1, seq2):

    return pose, movemap, dof_dict




def setup_dof_dict(seq1, seq2='')

    return dof_dict
