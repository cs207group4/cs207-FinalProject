class chemkin:

     """
    This class Chemkin computes the the reaction rates/ progress rates of the species 
    at each temperature of interest given species concentrations
    
    
    INPUTS
    =======
    Initialize with matrix of reactant, matrix of product and reaction_coeffs
    Using an XML file, the class calls InputParser and ReactionCoeffs to calculate the reaction rates

    RETURNS
    ========
   ###### 
   Returns the reaction rates/ progress rates of the species 
    at each temperature of interest given species concentrations
    
    EXAMPLES
    =========
     ## HOW CAN I PRINT AN EXAMPLE WITHOUT USING THE XML FILE? 
    ###>>> 

    
    
"""

    def __init__(self,nu_react,nu_prod,reaction_coeffs,species=None):
        self.nu_react = np.array(nu_react)
        self.nu_prod = np.array(nu_prod)
        if self.nu_prod.shape != self.nu_react.shape or len(reaction_coeffs) != self.nu_prod.shape[1]:
            raise ValueError("Dimensions not consistant!")
        self.rc_list = reaction_coeffs
        self.species = species

    @classmethod
    def from_xml(cls, filename):
        """
        calls Input Parser to parse xml file
        """
        input_ = InputParser(filename)
        rc_list = [ReactionCoeffs(**params) for params in input_.rate_coeff_params]
        return cls(input_.nu_react,input_.nu_prod,rc_list,input_.species)

    @classmethod
    def init_const_rc(cls,nu_react,nu_prod,rcs,species=None):
        rc_list = [ReactionCoeffs(type="Constant", k=rc_) for rc_ in rcs]
        return cls(nu_react,nu_prod,rc_list,species)

    def set_rc_params(self,**kwargs):
        for rc in self.rc_list:
            rc.set_params(**kwargs)

    def __repr__(self):
        class_name = type(self).__name__
        args = repr(self.nu_react.tolist()) + ',' + repr(self.nu_prod.tolist()) + ',' + repr(self.rc_list) + ',' + repr(self.species)
        return class_name + "(" + args +")"

    def __len__(self):
        '''return number of reactions'''
        return len(self.rc_list)

    def __str__(self):
        species_str = "chemkin: " + (str(self.species) if self.species else "")
        nu_react_str = "nu_react:\n" + str(self.nu_react)
        nu_prod_str = "nu_prod:\n" + str(self.nu_prod)
        rc_str = "reaction coefficients:\n[\n" + "\n".join([str(rc_) for rc_ in self.rc_list]) + "\n]"
        return "\n".join([species_str,nu_react_str,nu_prod_str,rc_str])

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