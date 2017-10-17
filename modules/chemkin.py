class chemkin:

    '''
    This class Chemkin computes the the reaction rates/ progress rates of the species
    at each temperature of interest given species concentrations


    INPUTS
    =======
    Initialize with matrix of reactant, matrix of product and reaction_coeffs
    Using an XML file, the class calls InputParser and ReactionCoeffs to calculate the reaction rates

    RETURNS
    ========
    After initialization, user could call:
     - set_rc_params(T=..., R=..., A=...): method to set params of reaction coeffs
     - reaction_rate(x): method to calculate reaction rate given concentration x
     - reaction_rate_T(x,T): method to calculate reaction rate given concentracion x and temprature T.
     - species: A variable containing species' names.
     - progress_rate(x): calculate progress rate given x...

    EXAMPLES
    =========
    >>> chem = chemkin.from_xml("rxns.xml")
    Finished reading xml input file
    >>> print(chem.species)
    ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    >>> chem.reaction_rate_T([1,1,1,1,1,1],1000)
    array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
            -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    '''

    def __init__(self,nu_react,nu_prod,reaction_coeffs,species=None, equations = None):
        self.nu_react = np.array(nu_react)
        self.nu_prod = np.array(nu_prod)
        if self.nu_prod.shape != self.nu_react.shape or len(reaction_coeffs) != self.nu_prod.shape[1]:
            raise ValueError("Dimensions not consistent!")
        self.rc_list = reaction_coeffs
        self.species = species
        self.equations = equations

    @classmethod
    def from_xml(cls, filename):
        """
        calls Input Parser to parse xml file, returns an initialized object
        """
        input_ = InputParser(filename)
        rc_list = [ReactionCoeffs(**params) for params in input_.rate_coeff_params]
        rxndata = input_.reactions
        equationlist = []
        for i, reaction in enumerate(rxndata):
            equationlist.append(rxndata[i]['equation'])
        return cls(input_.nu_react,input_.nu_prod,rc_list,input_.species, equationlist)

    @classmethod
    def init_const_rc(cls,nu_react,nu_prod,rcs,species=None):
        ''' construct an object with all constant or precalculated coeffs.
        '''
        rc_list = [ReactionCoeffs(type="Constant", k=rc_) for rc_ in rcs]
        return cls(nu_react,nu_prod,rc_list,species)

    def set_rc_params(self,**kwargs):
        ''' add new or change old parameters for all reaction coeffs.
        '''
        for rc in self.rc_list:
            rc.set_params(**kwargs)

    def __repr__(self):
        class_name = type(self).__name__
        args = repr(self.nu_react.tolist()) + ',' + repr(self.nu_prod.tolist()) + ',' + repr(self.rc_list) + ',' + repr(self.species)+ ',' + repr(self.equations)
        return class_name + "(" + args +")"

    def __len__(self):
        '''return number of reactions'''
        return len(self.rc_list)

    def __str__(self):
        if self.equations is not None:
            eqn_str = "chemical equations:\n[\n" + "\n".join([str(eq_) for eq_ in self.equations]) + "\n]"
        else:
            eqn_str = "chemical equations: not specified"
        if self.species is not None:

            species_str = "species: " + (str(self.species) if self.species else "")
        else:
            species_str = "species: not specified"
        nu_react_str = "nu_react:\n" + str(self.nu_react)
        nu_prod_str = "nu_prod:\n" + str(self.nu_prod)
        rc_str = "reaction coefficients:\n[\n" + "\n".join([str(rc_) for rc_ in self.rc_list]) + "\n]"
        return "\n".join([eqn_str, species_str,nu_react_str,nu_prod_str,rc_str])

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
        if len(x) != self.nu_prod.shape[0]:
            raise ValueError("ERROR: The concentration vector x must be of length i, where i is the number of species")
        #check that concentrations are all non-negative:
        if np.any(x)<0:
            raise ValueError("ERROR: All the species concentrations must be non-negative")
        #make the shape compatible with what NumPy needs for vectorized operations
        x = np.reshape(x, (len(x), 1))

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
        if len(x) != self.nu_prod.shape[0]:
            raise ValueError("ERROR: The concentration vector x must be of length i, where i is the number of species")
        #check that concentrations are all non-negative:
        if np.any(x)<0:
            raise ValueError("ERROR: All the species concentrations must be non-negative")
        #make the shape compatible with what NumPy needs for vectorized operations
        x = np.reshape(x, (len(x), 1))

        r = self.progress_rate(x)

        # Return an array...
        return np.sum(r * (self.nu_prod-self.nu_react), axis=1)

    def reaction_rate_T(self, x, T):
        '''
        A function to easily calculate reaction rate based on x and T.
        '''
        x = np.array(x)
        if len(x) != self.nu_prod.shape[0]:
            raise ValueError("ERROR: The concentration vector x must be of length i, where i is the number of species")
        #check that concentrations are all non-negative:
        if np.any(x)<0:
            raise ValueError("ERROR: All the species concentrations must be non-negative")
        #make the shape compatible with what NumPy needs for vectorized operations
        x = np.reshape(x, (len(x), 1))

        if T < 0:
            raise ValueError("ERROR: Temperature cannot be negative")
        self.set_rc_params(T=T)
        return self.reaction_rate(x)
