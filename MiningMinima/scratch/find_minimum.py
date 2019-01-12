def find_minimum(input_pose, scorefxn, movemap)
	'''Minimizes an input pose using indicated scorefunction and movemap'''
	
	minmover = rosetta.protocols.minimization_packing.MinMover(movemap, scorefxn,
			'lbfgs_armijo_nonmonotone', 1.0e-7 True)
	minmover.max_iter(1e7)
	minmover.min_options().use_nblist(True)
	minmover.min_options().nblist_auto_update(True)
	
	minmover.apply(input_pose)
	
	return input_pose