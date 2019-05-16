import numpy as np
import itertools
from hessian import numpy_hessian
from pyrosetta import *
from pyrosetta.rosetta import *
init()
def calc_hessian_at_min(pose_min, scorefxn, dof_dict):
    '''Take a minimized pose and computes the hessian at the base of the corresponding energy well.
	Returns the hessian matrix which may subsequently be diagonalized'''
	
    h = 0.1
    sig = 0.6
	
    h_rad = h*np.pi/180
    E_0 = scorefxn(pose_min)

    pose = Pose()
    pose.assign(pose_min)

    dofs = [pose_min.torsion(dof_dict[key]) for key in dof_dict]

    hessian_at_min = np.zeros((len(dofs), len(dofs)))

    n_pts = int(sig/h)
    if n_pts % 2 == 0: n_pts += 1
	
    ind = int(n_pts/2)

    tor_ranges = np.zeros((n_pts, len(dofs)))
    energy = np.zeros((n_pts, n_pts))

    for ii, dof in enumerate(dofs):

        tor_ranges[:, ii] = dof + np.linspace(-sig/2, sig/2, n_pts)

    for pair in list(itertools.combinations(dof_dict.keys(), 2)):

        foo = pair[0]
        bar = pair[1]

        x = tor_ranges[:, foo]
        y = tor_ranges[:, bar]

        for ii, xx in enumerate(x):

            pose.set_torsion(dof_dict[pair[0]], xx)

            for jj, yy in enumerate(y):

                pose.set_torsion(dof_dict[pair[1]], yy)
                energy[ii, jj] = scorefxn(pose)
                pose.set_torsion(dof_dict[pair[1]], dofs[bar])

            pose.set_torsion(dof_dict[pair[0]], dofs[foo])

        hess = numpy_hessian(energy, h_rad)

        #d2E_dx2 =  hess[0,0,:,:]
        #d2E_dy2 =  hess[1,1,:,:]
        d2E_dxdy = hess[0,1,:,:]

        #if hessian_at_min[foo, foo] == 0:

        #    hessian_at_min[foo, foo] = d2E_dx2[ind, ind]

        #if hessian_at_min[bar, bar] == 0:

        #    hessian_at_min[bar, bar] = d2E_dy2[ind, ind]

        hessian_at_min[foo, bar] = d2E_dxdy[ind, ind]
        hessian_at_min[bar, foo] = d2E_dxdy[ind, ind]
    
	pose.assign(pose_min)
	
    for key in dof_dict:
		
        energy = np.zeros((n_pts, 1))
        x = tor_ranges[:, key]
		
        for ii, xx in enumerate(x):
		
            pose.set_torsion(dof_dict[key], xx)
            energy[ii] = scorefxn(pose)
            pose.set_torsion(dof_dict[key], dofs[key])
			
	
        hessian_at_min[key, key] = (energy[ind+1] - 2*energy[ind] + energy[ind-1])/h_rad**2
		
    return hessian_at_min
