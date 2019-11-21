#!/bin/usr/env python

from recces.write_slurm_scripts import write_slurm_scripts
from recces.util import weight_evaluate

(temps, wts) = weight_evaluate('./prerun', './prerun_hist_scores.gz')

write_slurm_scripts(st_weights=wts,prerun=False)
