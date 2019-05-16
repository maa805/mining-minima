from pyrosetta import *
from pyrosetta.rosetta import *
from hessian import numpy_hessian
import numpy as np
import itertools
import pyrosetta.rosetta.protocols.rigid as rigid_moves

axis_dict = {3:numeric.xyzVector_double_t(1, 0, 0), 
4:numeric.xyzVector_double_t(0, 1, 0), 
5:numeric.xyzVector_double_t(0, 0, 1)}
origin = numeric.xyzVector_double_t(0,0,0)


def add_rb_dofs(pose, movemap, dof_dict):

	
	n_residues = pose.total_residue()
	start_res = pose.residue(1)
	end_res = pose.residue(n_residues)
	
	#ft = FoldTree( n_residues )
	ft = pose.fold_tree()
	ft.new_jump( 1, n_residues, n_residues/2 )
	
	ft.set_jump_atoms(1, start_res.atom_name(start_res.chi_atoms()[1][4]),
		end_res.atom_name(end_res.chi_atoms()[1][4]))
		
	pose.fold_tree(ft)
	
	movemap.set_jump(1, True)
	
	


# Need to figure out best way to deal with rb movers
# Rhiju accessed the rotation matrix and translation vector directly?
# Maybe a better strategy is to have a function to initialize all the relevant rb_movers
# Then call them as needed based on the dofs probed
# Ultimately, it may not matter as it essential comes down to a matter of bookkeeping

def initialize_rb_movers(pose, scorefxn, h=0.001):
	
	# Initialize three principle axes
	x_vec = numeric.xyzVector_double_t(1, 0, 0)
	y_vec = numeric.xyzVector_double_t(0, 1, 0)
	z_vec = numeric.xyzVector_double_t(0, 0, 1)
	
	# Set up translation movers
	step_x = rigid_moves.RigidBodyTransMover(pose, 1)	
	step_x.trans_axis(x_vec)
	step_x.step_size(h)
	
	step_y = rigid_moves.RigidBodyTransMover(pose, 1)
	step_y.trans_axis(y_vec)
	step_y.step_size(h)
	
	step_z = rigid_moves.RigidBodyTransMover(pose, 1)
	step_z.trans_axis(z_vec)
	step_z.step_size(h)
	
	# Set up rotation movers
	
	step_vx = rigid_moves.RigidBodyDeterministicSpinMover()
	step_vx.spin_axis(x_vec)
	#step_vx.rot_center(core.conformation.all_atom_center( pose.residue(2) ))
	step_vx.angle_magnitude(h*180./np.pi)
	
	step_vy = rigid_moves.RigidBodyDeterministicSpinMover()
	step_vy.spin_axis(y_vec)
	#step_vy.rot_center(core.conformation.all_atom_center( pose.residue(2) ))
	step_vy.angle_magnitude(h*180./np.pi)
	
	step_vz = rigid_moves.RigidBodyDeterministicSpinMover()
	step_vz.spin_axis(z_vec)
	#step_vz.rot_center(core.conformation.all_atom_center( pose.residue(2) ))
	step_vz.angle_magnitude(h*180./np.pi)
	
	return (step_x, step_y, step_z, step_vx, step_vy, step_vz)
	
	
def calc_rb_hessian(pose, scorefxn, rb_movers, h=0.001, w = 0.006):
	
	hessian = np.zeros( (6,6) )
	
	n_pts = int(w/h)
	if n_pts % 2 == 0: n_pts += 1
	ind = int(n_pts/2)	
	
	minimum = Pose()
	minimum.assign(pose)
	 
	for pair in list(itertools.combinations(range(6), 2)):
	
		dof_1 = pair[0]
		dof_2 = pair[1]
		energy = rb_energy(pose, scorefxn, rb_movers, dof_1, dof_2, h)
		pose.assign(minimum)
	 
		hess = numpy_hessian(energy, h=h)
    
		#d2E_dx2 = hess[0,0,:,:]
		#d2E_dy2 = hess[1,1,:,:]
		d2E_dxdy = hess[0,1,:,:]
    
		#if hessian[dof_1, dof_1] == 0: hessian[dof_1, dof_1] = d2E_dx2[ind, ind]
		#if hessian[dof_2, dof_2] == 0: hessian[dof_2, dof_2] = d2E_dy2[ind, ind]
		hessian[dof_1, dof_2] = d2E_dxdy[ind, ind]
		hessian[dof_2, dof_1] = d2E_dxdy[ind, ind]
		
	for dof in range(6):
		
		energy = rb_energy(pose, scorefxn, rb_movers, dof, dof)
		pose.assign(minimum)
		
		hess = (energy[ind-1] - 2*energy[ind] + energy[ind+1])/h**2
		
		hessian[dof, dof] = hess
    
	pose.assign(minimum)
	return hessian
	
	
def rb_energy(pose, scorefxn, rb_movers, dof_1, dof_2, h=0.001, w=0.006):
    
	n_pts = int(w/h)
	if n_pts % 2 == 0: n_pts += 1
	
	dof_1_mover = rb_movers[dof_1]

	# Set up start and reset movers for current dofs	
	dof_1_start_mover = set_up_start_mover(dof_1, w, h, pose)
	
	dof_1_start_mover.apply(pose)

	if dof_1 == dof_2:
	
		energy = np.zeros( (n_pts, 1) )
		
		dof_1_start_mover.apply(pose)
		
		energy[0] = scorefxn(pose)
		
		for ddof in range(1, n_pts):
		
			dof_1_mover.apply(pose)
			energy[ddof] = scorefxn(pose)
			
	else:
	
		dof_2_mover = rb_movers[dof_2]
		dof_2_start_mover = set_up_start_mover(dof_2, w, h, pose)
		dof_2_reset_mover = set_up_reset_mover(dof_2, w, h, pose)
		
		energy = np.zeros( (n_pts, n_pts) )
		
		dof_2_start_mover.apply(pose)
			
		energy[0,0] = scorefxn(pose)
		
		for ddof_2 in range(1, n_pts):
			
			dof_2_mover.apply(pose)
			energy[0, ddof_2] = scorefxn(pose)
			
		dof_2_reset_mover.apply(pose)
		
		for ddof_1 in range(1, n_pts):
			
			dof_1_mover.apply(pose)
			
			for ddof_2 in range(n_pts):
				
				dof_2_mover.apply(pose)
				energy[ddof_1, ddof_2] = scorefxn(pose)
			
			dof_2_reset_mover.apply(pose)
		
	return energy


def set_up_start_mover(dof, w, h, pose):

	if dof in [0, 1, 2]:
	
		start_mover = rigid_moves.RigidBodyTransMover()
		start_mover.trans_axis( axis_dict[dof+3] )
		start_mover.step_size( -w/2 )
			
	else:
	
		start_mover = rigid_moves.RigidBodyDeterministicSpinMover()
		#start_mover.rot_center( core.conformation.all_atom_center( pose.residue(2) ) )
		start_mover.spin_axis( axis_dict[dof] )
		start_mover.angle_magnitude(-w/2./np.pi*180.) 
	
	return start_mover
	

def set_up_reset_mover(dof, w, h, pose):

	if dof in [0, 1, 2]:
	
		reset_mover = rigid_moves.RigidBodyTransMover()
		reset_mover.trans_axis( axis_dict[dof+3] )
		reset_mover.step_size( -(w+h) )
			
	else:
	
		reset_mover = rigid_moves.RigidBodyDeterministicSpinMover()
		#reset_mover.rot_center( core.conformation.all_atom_center( pose.residue(2) ) )
		reset_mover.spin_axis( axis_dict[dof] )
		reset_mover.angle_magnitude(-(w+h)/np.pi*180.) 
		
	return reset_mover
	
def rb_mode_scan(min_pose, rb_movers, mode, scorefxn, h = 0.02, w = 6.0, kt = 1.0):
	'''
	Given an input pose (assumed to be a local minimum) and a correspondig mode, returns the partition
	function evaluated along the mode.'''

	# Initialize stuff 
	min_energy = scorefxn(min_pose)
	energy = min_energy
	
	pose = Pose()
	pose.assign(min_pose)
	obs = protocols.moves.AddPyMOLObserver(pose, keep_history=True)
	
	for dof in range(6):

		if dof in [0, 1, 2]: rb_movers[dof].step_size(h*mode[dof])	
		else: rb_movers[dof].angle_magnitude(h*mode[dof]/np.pi*180.)
			
	n_pts = int(w/h + 1)
	result  = [scorefxn(pose)]
	
	# Scan along mode at indicated granularity 
	for ii in range(n_pts/2):
	
		for dof in range(6): rb_movers[dof].apply(pose)
		result.append( scorefxn(pose) )
	
	pose.assign(min_pose)
	
	for dof in range(6):

		if dof in [0, 1, 2]: rb_movers[dof].step_size(-h*mode[dof])	
		else: rb_movers[dof].angle_magnitude(-h*mode[dof]/np.pi*180.)	

	for ii in range(n_pts/2):
	
		for dof in range(6): rb_movers[dof].apply(pose)
		result.insert(0, scorefxn(pose) )
	
	
	# Apply trapezoidal rule to get anharmonic free energy
	partition_function = np.trapz( np.exp( -( np.array(result) - scorefxn(min_pose) )/kt ), dx=h)
	
	return partition_function, result

	
	
	
