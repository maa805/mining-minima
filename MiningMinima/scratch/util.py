import numpy as np

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def within_360(dof):
    
    diff = dof/360.0
    diff = np.round(diff)
    
    return dof - 360.0*diff
	
def is_repeat_minimum(new_dofs, dof_lists):
    
    #dofs_rounded = np.around(new_dofs, decimals=1)
    
    if len(dof_lists) == 0:
        
        return False 
    
    for dof_list in dof_lists:
        
        #dof_list_rounded = np.around(dof_list, decimals=1)

        dof_diff = within_360(new_dofs) - within_360(dof_list)
        
        if np.all(abs(within_360(dof_diff)) < 20.0): return True
        else: continue
    
    return False