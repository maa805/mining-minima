from pyrosetta import *
from pyrosetta.rosetta import *
from hessian import numpy_hessian
import numpy as np
import itertools
import pyrosetta.rosetta.protocols.rigid as rigid_moves

axis_dict = {3:numeric.xyzVector_double_t(1, 0, 0), 
4:numeric.xyzVector_double_t(0, 1, 0), 
5:numeric.xyzVector_double_t(0, 0, 1)}


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
	step_vx.rot_center(core.pose.center_of_mass(pose, 1, 2))
	
	step_vy = rigid_moves.RigidBodyDeterministicSpinMover()
	step_vy.spin_axis(y_vec)
	step_vy.rot_center(core.pose.center_of_mass(pose, 1, 2))
	
	step_vz = rigid_moves.RigidBodyDeterministicSpinMover()
	step_vz.spin_axis(z_vec)
	step_vz.rot_center(core.pose.center_of_mass(pose, 1, 2))
	
	return (step_x, step_y, step_z, step_vx, step_vy, step_vz)
	
	
def calc_rb_hessian(pose, scorefxn, rb_movers, h=0.001):
	
	hessian = np.zeros( (6,6) )
	
	w = 0.004
	n_pts = int(w/h + 1)
	ind = int(n_pts/2)
	
	minimum = Pose()
	minimum.assign(pose)
	 
	for pair in list(itertools.combinations(range(6), 2)):
	
		i = pair[0]
		j = pair[1]
		energy = rb_energy(pose, scorefxn, rb_movers, i, j, h)
		pose.assign(minimum)
	 
		hess = numpy_hessian(energy, h=h)
    
		d2E_dx2 = hess[0,0,:,:]
		d2E_dy2 = hess[1,1,:,:]
		d2E_dxdy = hess[0,1,:,:]
    
		if hessian[i, i] == 0: hessian[i, i] = d2E_dx2[ind, ind]
		if hessian[j, j] == 0: hessian[j, j] = d2E_dy2[ind, ind]
		hessian[i, j] = d2E_dxdy[ind, ind]
		hessian[j, i] = d2E_dxdy[ind, ind]
    
	pose.assign(minimum)
	return hessian
	
	
def rb_energy(pose, scorefxn, rb_movers, i, j, h=0.001):
    
	w = 0.004
	n_pts = int(w/h + 1)
	energy = np.zeros( (n_pts,n_pts) )
    
	i_mover = rb_movers[i]
	j_mover = rb_movers[j]
	
	if i in [0, 1, 2]: 
		
		i_start_mover = rigid_moves.RigidBodyTransMover()
		i_start_mover.trans_axis( i_mover.trans_axis() )
		i_start_mover.step_size(-w/2)
		
	else: 
		
		i_start_mover = rigid_moves.RigidBodyDeterministicSpinMover()
		i_start_mover.spin_axis( axis_dict[i] )
		i_start_mover.rot_center( core.pose.center_of_mass(pose, 1, 2) )
		i_start_mover.angle_magnitude(-w/2./np.pi*180.)
		
	if j in [0, 1, 2]: 
		
		j_start_mover = rigid_moves.RigidBodyTransMover()
		j_start_mover.trans_axis( j_mover.trans_axis() )
		j_start_mover.step_size(-w/2)
		
		j_reset_mover = rigid_moves.RigidBodyTransMover()
		j_reset_mover.trans_axis( j_mover.trans_axis() )
		j_reset_mover.step_size( -(w+h) )
		
	else: 
		
		j_start_mover = rigid_moves.RigidBodyDeterministicSpinMover()
		j_start_mover.spin_axis( axis_dict[j] )
		j_start_mover.rot_center( core.pose.center_of_mass(pose, 1, 2) )
		j_start_mover.angle_magnitude(-w/2./np.pi*180.)
		
		j_reset_mover = rigid_moves.RigidBodyDeterministicSpinMover()
		j_reset_mover.spin_axis( axis_dict[j] )
		j_reset_mover.rot_center( core.pose.center_of_mass(pose, 1, 2) )
		j_reset_mover.angle_magnitude(-(w+h)/np.pi*180.)
		    
	i_start_mover.apply(pose)
	j_start_mover.apply(pose)
	
	energy[0,0] = scorefxn(pose)

	for dy in range(1, n_pts):
        
		j_mover.apply(pose)
		energy[0, dy] = scorefxn(pose)
        
	j_reset_mover.apply(pose)
	
	for dx in range(1, n_pts):
        
		i_mover.apply(pose)
        
		for dy in range(n_pts):
            
			j_mover.apply(pose)
			energy[dx,dy] = scorefxn(pose)
		
		j_reset_mover.apply(pose)
		
	return energy


	
	
	
	
