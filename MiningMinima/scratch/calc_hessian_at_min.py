import numpy as np
import itertools
from hessian import numpy_hessian
from pyrosetta import *
from pyrosetta.rosetta import *
init()
def calc_hessian_at_min(pose_min, scorefxn, dof_dict):
    '''Take a minimized pose and computes the hessian at the base of the corresponding energy well.
	Returns the hessian matrix which may subsequently be diagonalized'''
    print dof_dict
    h = 0.1
    h_rad = h*np.pi/180
    E_0 = scorefxn(pose_min)

    pose = Pose()
    pose.assign(pose_min)

    dofs = [pose_min.torsion(dof_dict[key]) for key in dof_dict]

    hessian_at_min = np.zeros((len(dofs), len(dofs)))

    n_pts = int(1/h + 1)
    ind = int(n_pts/2)

    tor_ranges = np.zeros((n_pts, len(dofs)))
    energy = np.zeros((n_pts, n_pts))

    for ii, dof in enumerate(dofs):

        tor_ranges[:, ii] = dof + np.arange(-0.5, 0.5+h, h)

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

        d2E_dx2 =  hess[0,0,:,:]
        d2E_dy2 =  hess[1,1,:,:]
        d2E_dxdy = hess[0,1,:,:]

        if hessian_at_min[foo, foo] == 0:

            hessian_at_min[foo, foo] = d2E_dx2[ind, ind]

        if hessian_at_min[bar, bar] == 0:

            hessian_at_min[bar, bar] = d2E_dy2[ind, ind]

        hessian_at_min[foo, bar] = d2E_dxdy[ind, ind]
        hessian_at_min[bar, foo] = d2E_dxdy[ind, ind]
        
    return hessian_at_min
