from pyrosetta import *
from pyrosetta.rosetta import *
import numpy as np

def set_up_multifunc(pose, movemap, scorefxn):
    # set up minimizer map
    start_score = scorefxn(pose)
    min_map = core.optimization.MinimizerMap()
    min_map.setup(pose, movemap)
    pose.energies().set_use_nblist(pose, min_map.domain_map(), True)
    scorefxn.setup_for_minimizing(pose, min_map)
    multifunc = core.optimization.AtomTreeMultifunc(pose, min_map, scorefxn)

    # values to return 
    return min_map, multifunc

def find_minimum(min_pose, min_map, multifunc, scorefxn, min_options=False):
    # use default minimization options unless user supplies own
    minimizer_options = core.optimization.MinimizerOptions(
        'lbfgs_armijo_nonmonotone', 1.5e-9, True, False, False)
    minimizer_options.max_iter(10000)
    
    # add support for user supplied min options
    
    start_score = scorefxn(min_pose)
    min_pose.energies().set_use_nblist(min_pose, min_map.domain_map(), True)
    scorefxn.setup_for_minimizing(min_pose, min_map)
    scorefxn.setup_for_derivatives(min_pose)
    
    minimizer_options.use_nblist(True)
    minimizer_options.nblist_auto_update(True) # for some reason off by default 

    # initialize min_dofs and gradient containers
    min_dofs = Vector1([0.0]*min_map.nangles())

    # get dofs from pose
    min_map.copy_dofs_from_pose(min_pose, min_dofs)

    # set up minimizer and run
    minimizer = core.optimization.Minimizer(multifunc, minimizer_options)
    minimizer.run(min_dofs)

    # copy new min dofs to min pose
    min_map.copy_dofs_to_pose(min_pose, min_dofs)
    return np.array(min_dofs)