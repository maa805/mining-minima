import sys
import os

def write_slurm_scripts(st_weights='', temps=[0.5, 0.8, 1.0, 1.4, 1.8, 3.0, 7.0, 30.0], prerun=True):
	dir_name = os.path.basename( os.getcwd() )
	
	seqs = dir_name.split('_')
	seq1 = seqs[0]
	if len(seqs) == 2: seq2 = seqs[-1]
	else: seq2 = ''

	if prerun:
	
		os.chdir('prerun/')
		f = open('prerun.sh', 'w+')
		
	else:
		
		os.chdir('ST/')
		f = open('ST.sh', 'w+')
		
	f.write('#!/bin/bash\n')
	f.write('#SBATCH --job-name=' + seq1 + seq2 + '\n')
	f.write('#SBATCH --mem-per-cpu=1G\n')
	f.write('#SBATCH --cpus-per-task=1\n')
	f.write('#SBATCH --partition=biochem\n')
	f.write('#SBATCH --mail-type=FAIL\n')
	
	base_cmd_line = 'srun --ntasks=1 $ROSETTA/main/source/bin/recces.linuxgccrelease -score:weights stepwise/rna/turner_no_zeros '
	
	base_cmd_line += ' -seq1 ' + seq1
	
	if seq2: base_cmd_line += ' -seq2 ' + seq2
	
	if prerun:
	
		f.write('#SBATCH --time=06:00:00\n')
		f.write('#SBATCH --ntasks=8\n')
		f.write('\n')
		#f.write('ml biology rosetta\n')
		
		base_cmd_line += ' -n_cycle 300000 -out_prefix prerun'
		
		for temp in temps: f.write(base_cmd_line + ' -temps ' + str(temp) + ' &\n')
			
	else:
	
		temps_str = ''
		st_weights_str = ''
		
		for temp in temps: temps_str += str(temp) + ' '
		for weight in st_weights: st_weights_str += str(weight) + ' '
		
		base_cmd_line += ' -save_score_terms '
		
		f.write('#SBATCH --time=12:00:00\n')
		f.write('#SBATCH --ntasks=2\n')
		f.write('\n')
		#f.write('ml biology rosetta \n')
		
		f.write(base_cmd_line + ' -temps ' + temps_str + ' -st_weights ' + st_weights_str + 
				' -n_cycle 9000000 -out_prefix ST -dump_pdb -dump_freq 9000 &\n')
		f.write(base_cmd_line + ' -temps -1 -n_cycle 300000 -out_prefix kT_inf &\n')
		
	f.write('wait\n')
	f.close()


