import os
import argparse

parser = argparse.ArgumentParser(description='Make default directories for RECCES-style simulations.')
parser.add_argument('--sim_file', metavar='infile', type=str)
parser.add_argument('--simulations', metavar='sim', type=str, nargs ='+')

args = parser.parse_args()

if args.sim_file:
	sim_list = [line.split() for line in open(args.sim_file)]
	sim_list = [sim_list[ii-1][0] for ii in range(len(sim_list))]
else:
	sim_list = args.simulations

for sim in sim_list:
	
	try: os.mkdir(sim)
	except OSError: pass

	os.chdir(sim)
	
	try: os.mkdir('prerun/')
	except OSError: pass

	try: os.mkdir('ST/')
	except OSError: pass
	
	os.chdir('./..')
