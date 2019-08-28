import numpy as np
import scipy
import scipy.special
import itertools

from pyrosetta import *
from pyrosetta.rosetta import *
init()

from pose_setup import add_bb_suite
from pose_setup import add_chi_dofs

KT_IN_KCAL = 0.6163

class MiningMinima:	
    def __init__(self, seq1='', seq2='', infile='', scorefxn = 'stepwise/rna/turner', kt = 1.0):
        self.kt = kt
		
        self.harmonic_free_energy = 0.0
        self.anharmonic_free_energy = 0.0
		
        self.input_pose = Pose()
        if not infile: 
			
            self.seq1 = seq1
            self.seq2 = seq2
            self.pose_setup_turner()
			
        else: 
            self.infile = infile
            self.pose_setup_from_file()
		
        self.n_dofs = len(self.dof_dict)
		
        self.scorefxn = core.scoring.ScoreFunctionFactory.create_score_function(scorefxn)
		
        # set up the multifunc and associated acts
        self.min_pose = Pose()
        self.min_pose.assign(self.input_pose)
        self.set_up_multifunc(self.min_pose)
        
        # find the min_dofs
        self.min_dofs = self.find_minimum()
        self.min_energy = self.scorefxn(self.min_pose)
        
        # Calculate hessian at base of minimum and normal modes
        self.hessian = hessian_at_min(self.min_dofs, self.multifunc)
        self.eigenvalues, self.modes = np.linalg.eigh(self.hessian)
		
		# Function for density of states
        #self.dos = lambda E: (2.0*np.pi)**(self.n_dofs/2)*(E-self.min_energy)**(self.n_dofs/2 -1 )/scipy.special.gamma(self.n_dofs/2)/np.sqrt(np.linalg.det(self.hessian))*np.heaviside(E-self.min_energy, 0.5)
		#self.dos = lambda E: 0.5/self.n_dofs*(2.0*np.pi*(E-self.min_energy))**(self.n_dofs/2)/scipy.special.gamma(self.n_dofs/2)/np.sqrt(np.linalg.det(self.hessian))*np.heaviside(E-self.min_energy, 0.5)
		
		# Calculate free energy
        #self.calc_harmonic_free_energy()
        #self.calc_anharmonic_free_energy()
        
    def pose_setup_turner(self):
        n_residues = len(self.seq1)
        if self.seq2: n_residues += len(self.seq2)
        
        dof_dict = {}
        movemap = MoveMap()
		
        pose = protocols.recces.pose_setup_turner(self.seq1, self.seq2)
		
        idx = 0
        suite = 2
		
        for i in range(n_residues - 1):
            if pose.fold_tree().is_cutpoint(i+1):
                suite += 1
                continue 
			
            movemap, dof_dict = add_bb_suite(suite, idx, movemap, dof_dict)				
            idx += 1
            suite += 1			
		
        movemap, dof_dict = add_chi_dofs(n_residues, idx, movemap, dof_dict)

        self.input_pose.assign(pose)
        self.movemap = movemap
        self.dof_dict = dof_dict
		
	
    def pose_setup_from_file(self): 
        pose = pose_from_file(self.infile)
        n_residues = len(pose.sequence())
        
        ft = FoldTree(n_residues)
        ft.new_jump(1, n_residues, n_residues/2)
        # use usual chi torsion atoms to set jump
        ft.set_jump_atoms(1, pose.residue(1).atom_name(pose.residue(1).chi_atoms(1)[4]),
            pose.residue(n_residues).atom_name(pose.residue(n_residues).chi_atoms(1)[4]))
		
        pose.fold_tree(ft)
	
        dof_dict = {}
        movemap = MoveMap()
		
        idx = 0
        suite = 2
		
        for i in range(n_residues - 1):
		
            if pose.fold_tree().is_cutpoint(i+1):
                suite += 1
                continue
			
            movemap, dof_dict = add_bb_suite(suite, idx, movemap, dof_dict)
            idx += 1
            suite += 1
			
        movemap, dof_dict = add_chi_dofs(n_residues, idx, movemap, dof_dict)
		
        self.input_pose.assign(pose)
        self.movemap = movemap
        self.dof_dict = dof_dict
		
        return
		
	#def calc_hessian_at_min(self, self.min_pose, self.scorefxn, self.dof_dict)
	
	
	#def mode_scan(self, self.min_pose, self.scorefxn, self.dof_dict) 
	
    def set_up_multifunc(self, pose):
        # set up minimizer map
        min_map = core.optimization.MinimizerMap()
        min_map.setup(pose,self.movemap)
        
        start_score = self.scorefxn(pose)
        self.input_pose.energies().set_use_nblist(pose,min_map.domain_map(),True)
        multifunc = core.optimization.AtomTreeMultifunc(pose,min_map,self.scorefxn)
        self.scorefxn.setup_for_minimizing(pose,min_map)
        self.scorefxn.setup_for_derivatives(pose)
   
        self.min_map = min_map
        self.multifunc = multifunc

    def find_minimum(self):
        # set min options
        min_options = core.optimization.MinimizerOptions(
            'lbfgs_armijo_nonmonotone', 1e-15, True, False, False)
        min_options.use_nblist(True)
        min_options.nblist_auto_update(True) # for some reason off by default 
        
        # initialize min_dofs and gradient containers
        min_dofs = Vector1([0.0]*self.n_dofs)
            
        # get dofs from pose
        self.min_map.copy_dofs_from_pose(self.min_pose, min_dofs)
            
        # set up minimizer and run
        minimizer = core.optimization.Minimizer(self.multifunc, min_options)
        for _ in range(1): minimizer.run(min_dofs)
        
        # copy new min dofs to min pose
        self.min_map.copy_dofs_to_pose(self.min_pose, min_dofs)
        return np.array(min_dofs)
    
    
# Older version of find_minimum using minmover  
#    def find_minimum(self):
#        '''Minimizes an input pose using indicated scorefunction and movemap''' 
#        minmover = rosetta.protocols.minimization_packing.MinMover(
#            self.movemap, self.scorefxn,'lbfgs_armijo_nonmonotone', 1.0e-7, True)
#        minmover.max_iter(1000000)
#        minmover.min_options().use_nblist(True)
#        minmover.min_options().nblist_auto_update(True)
#		
#        minmover.apply(self.min_pose)
#        self.min_map.setup(self.min_pose,self.movemap)
#        self.min_pose.energies().set_use_nblist(self.min_pose,self.min_map.domain_map(),True)
#        self.scorefxn.setup_for_minimizing(self.min_pose,self.min_map)
#        self.scorefxn.setup_for_derivatives(self.min_pose)
		
    def calc_anharmonic_free_energy(self):
        total_partition = 1.0
	
        for i in range(self.n_dofs): 
            mode_partition, _ = mode_scan(self.min_pose, self.modes[:,i], self.scorefxn, self.dof_dict)
            total_partition *= mode_partition
        
        self.anharmonic_free_energy = self.min_energy - self.kt*np.log(total_partition)

	def calc_harmonic_free_energy(self):
		self.harmonic_free_energy = (self.min_energy - 0.5*self.kt*self.n_dofs*np.log(2*np.pi*self.kt) + 0.5*self.kt*np.log(np.prod(self.eigenvalues)))
	
	def harmonic_ensemble(self, n_struct=200):
		'''
		Parameters
		--------------------
		n_struct:	Number of ensemble members to generate
		'''
		ensemble = np.zeros((self.n_dofs, n_struct))
		
		mu = [self.min_pose.torsion(self.dof_dict[key]) for key in self.dof_dict]
		cov = np.matmul( np.matmul( self.modes, np.diag(1/self.eigenvalues)), self.modes.T)*(180./np.pi)**2	
		
		ensemble = np.random.multivariate_normal(mu, cov, size=(n_struct))
		
		return ensemble
        
def hessian_at_min(min_dofs, multifunc, h=1e-3):
    min_dofs = np.array(min_dofs) # take advantage of numpy indexing
    n_dofs = len(min_dofs)
    hessian = np.zeros((n_dofs,n_dofs))
    plus = Vector1([0.0]*n_dofs)
    minus = Vector1([0.0]*n_dofs)
    
    for ii in range(n_dofs):
        new_dofs = min_dofs[:]
        new_dofs[ii] += h*180./np.pi # dofs in degrees
        multifunc.dfunc(Vector1(list(new_dofs)), plus)
        new_dofs[ii] -= 2.*h*180./np.pi
        multifunc.dfunc(Vector1(list(new_dofs)), minus)
        
        #central difference scheme
        row = (np.array(plus) - np.array(minus))/2./h
        hessian[ii] = row * 180./np.pi
    return 0.5*(hessian + hessian.T) # enforce symmetry
    
    
def mode_scan(min_dofs, multifunc, mode, limit=np.pi/3, dx=0.005):
    # convert to np array for useful indexing
    min_dofs = np.array(min_dofs)
    displacement_array = np.linspace(-limit, limit, int(2*limit/dx)+1)
    result = np.zeros_like(displacement_array)
    
    new_dofs = np.array(min_dofs)
    for ii, displacement in enumerate(displacement_array):
        new_dofs = min_dofs[:] + displacement*mode*180./np.pi
        result[ii] = multifunc(Vector1(list(new_dofs)))
 
    return result
 
 
def compute_total_partition(min_dofs, multifunc, modes, eigenvalues, limit=np.pi/3, dx=0.005):
    from scipy.special import erf
    
    total_log_partition = 0.
    total_log_harmonic = 0.
    scans = tuple()
    xx = np.linspace(-limit, limit, int(2*limit/dx)+1)
    min_energy = multifunc(Vector1(list(min_dofs)))

    for ii, mode in enumerate(modes.T):      # columns of array are eigenvectors
        result = mode_scan(min_dofs, multifunc, mode, limit=limit, dx=dx)
        result -= min_energy                 # energy is relative to base of well
        scans = scans + (result,)
        total_log_partition += np.log(np.trapz(np.exp(-result),dx=dx))
        total_log_harmonic += np.log(np.sqrt(2*np.pi/eigenvalues[ii])*erf(2*limit/np.sqrt(2/eigenvalues[ii])))

    return total_log_partition, total_log_harmonic, scans
    
def array_to_vector1(array):
    vector_from_array = Vector1(list(array))
    return vector_from_array