import numpy as np

class reaction_rates:
    '''
    Initialize with matrix of reactant, matrix of product and reaction_coeffs
    '''

    def __init__(self,nu_react,nu_prod,reaction_coeffs):
        self.nu_react = np.array(nu_react)
        self.nu_prod = np.array(nu_prod)
        if nu_prod.shape != nu_react.shape or len(reaction_coeffs) != nu_prod.shape[1]:
            raise ValueError("Dimensions not consistant!")
        self.rc_list = reaction_coeffs

    def set_rc_params(self,**kwargs):
        for rc in self.rc_list:
            rc.set_params(**kwargs)


    def progress_rate(self,x):
        ''' 
        return progress rate for reactions of form:
        v_11 A + v_21 B -> v_31 C
        v_12 A + v_32 C -> v_22 B + v_32 C
        
        or a more general form
        
        INPUTS
        =======
        x: A i*1 vector specifying concentration of each specie
        
        RETURNS
        ========
        R: Progress rates of j reactions (np.array)
        
        '''
        
        x = np.array(x)
        
        if x.shape[1] != 1 or x.shape[0] != self.nu_prod.shape[0]:
            raise ValueError("Must satisfy: x -> i*1 matrix, i is the number of species")
            
        # Return an array
        return np.array([rc.kval() for rc in self.rc_list]).astype(float) * np.product(x ** self.nu_react, axis = 0)

    def reaction_rate(self,x):
        ''' 
        return reaction rate for each species in the system:
        v_11 A + v_21 B -> v_31 C
        v_32 C -> v_12 A + v_22 B
        
        or a more general form
        
        INPUTS
        =======
        x: A i*1 vector specifying concentration of each specie
        
        RETURNS
        ========
        R: reaction rates of species (np.array)
        
        '''
        x = np.array(x)
        
        if x.shape[1] != 1 or x.shape[0] != self.nu_prod.shape[0]:
            raise ValueError("Must satisfy: x -> i*1 matrix, i is the number of species")
        
        r = self.progress_rate(x)
            
        # Return an array...
        return np.sum(r * (self.nu_prod-self.nu_react), axis=1)