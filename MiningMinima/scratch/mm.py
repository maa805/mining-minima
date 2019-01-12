from pyrosetta import *
from pyrosetta.rosetta import *
import argparse	

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates pose either from file or sequence')
    parser.add_argument('-in_file', help='Input pdb file to be processed')
    parser.add_argument('seq1', help='Sequence of first strand')
    parser.add_argument('seq2', help='Sequence of second strand')

    args = parser.parse_args()

    if args.infile:

         pose, movemap, dof_dict = pose_setup_from_file(in_file)

    if seq1:

         pose, movemap, dof_dict = pose_setup_turner(seq1, seq2)


