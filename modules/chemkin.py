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
     - reaction_rate(x): method to calculate reaction rate given concentration x (assumes temperature is already set)
     - reaction_rate_T(x,T): method to calculate reaction rate given concentration x and temperature T.
     - species: A sequence (list, array, tuple) containing sthe names of the species
     - progress_rate(x): calculate progress rate given x (assumes temperature is already set)

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

    def __init__(self,nu_react,nu_prod,reaction_coeffs,rxn_types, reversible, species = None, equations = None):
        """
        INPUT
        =====
        nu_prod: Array of integers, required
                 N X M array of stoichiometric coefficients for products (N species, M reactions)
        nu_react: Array of integers, required
                 N x M array of stoichiometric coefficients for reactants
        reaction_coeffs: Array, required
                 Length M list of reaction rate coefficients in the form of representations of ReactionCoeff objects
        rxn_types: Array of strings, required
                 List or array of length M of reaction types
        reversible: Array, required
                 List or array of length M denoting whether the reaction is reversible
        species: Array of strings, optional
                 List or array of length N providing the species
        equations: Array of strings, optional
                 List or array of length M providing the reaction equations
        """
        self.nu_react = np.array(nu_react)
        self.nu_prod = np.array(nu_prod)
        if self.nu_prod.shape != self.nu_react.shape or len(reaction_coeffs) != self.nu_prod.shape[1]:
            raise ValueError("Dimensions not consistent!")
        self.rc_list = reaction_coeffs
        self.species = species
        self.equations = equations
        self.rxn_types = rxn_types
        self.reversible = reversible

    @classmethod
    def from_xml(cls, filename):
        """
        Alternate constructor starting with an xml input file

        INPUT
        =====
        filename: string, required
                  Name of input xml file

        RETURNS
        =======
        initialized chemkin object
        """
        input_ = InputParser(filename)
        rc_list = [ReactionCoeffs(**params) for params in input_.rate_coeff_params]
        rxndata = input_.reactions
        equationlist = []
        rxn_types = []
        reversible = []
        for i, reaction in enumerate(rxndata):
            equationlist.append(reaction['equation'])
            rxn_types.append(reaction['type'])
            reversible.append(reaction['reversible'])
        return cls(input_.nu_react,input_.nu_prod,rc_list,rxn_types, reversible, species = input_.species, equations = equationlist)

    @classmethod
    def init_const_rc(cls,nu_react,nu_prod,rcs,rxn_types,reversible,species=None,equations=None):
        ''' construct an object with all constant or precalculated coeffs.
        '''
        rc_list = [ReactionCoeffs(type="Constant", k=rc_) for rc_ in rcs]
        return cls(nu_react,nu_prod,rc_list,rxn_types, reversible, species,equations)

    def set_rc_params(self,**kwargs):
        ''' add new or change old parameters for all reaction coeffs.
        '''
        for rc in self.rc_list:
            rc.set_params(**kwargs)

    def __repr__(self):
        class_name = type(self).__name__
        args = repr(self.nu_react.tolist()) + ',' + repr(self.nu_prod.tolist()) + ',' + repr(self.rc_list) + ',' + repr(self.rxn_types)+ ',' + repr(self.reversible)+',' + repr(self.species)+ ',' + repr(self.equations)
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
        rxn_str = "reaction types: " + str(self.rxn_types)
        reversible_str = "reversible: " + str(self.reversible)
        return "\n".join([eqn_str, species_str,nu_react_str,nu_prod_str,rc_str, rxn_str, reversible_str])

    def progress_rate(self,x):
        '''
        Return progress rate for a system of M reactions involving N species


        INPUTS
        =======
        x: array or list, required
           A length-N vector specifying concentration of each specie

        RETURNS
        ========
        R: Array of reaction rates of species (length N)

        '''

        x = np.array(x)
        if len(x) != self.nu_prod.shape[0]:
            raise ValueError("ERROR: The concentration vector x must be of length N, where N is the number of species")
        #check that concentrations are all non-negative:
        if np.any(x<0):
            raise ValueError("ERROR: All the species concentrations must be non-negative")
        #make the shape compatible with what NumPy needs for vectorized operations
        x = np.reshape(x, (len(x), 1))

        return np.array([rc.kval() for rc in self.rc_list]).astype(float) * np.product(x ** self.nu_react, axis = 0)

    def reaction_rate(self,x):
        '''
        Return reaction rates for a system of M reactions involving N species

        INPUTS
        =======
        x: array or list, required
           A length-N vector specifying concentration of each specie

        RETURNS
        ========
        R: Array of reaction rates of species (length N)

        '''

        r = self.progress_rate(x)

        # Return an array...
        return np.sum(r * (self.nu_prod-self.nu_react), axis=1)

    def reaction_rate_T(self, x, T):
        """
        A wrapper function to calculate reaction rate based on user-defined x and T

        INPUT
        ====
        x: Array or list, required
           Vector of concentration of each species of length N, where N is the number of species
        T: float, required
           Temperature, in Kelvin

        RETURNS
        ======
        Array of length N containing reaction rates for each species

        """

        if T <= 0:
            raise ValueError("ERROR: Temperature must be positive")
        self.set_rc_params(T=T)
        if np.any(np.array(self.rxn_types)!='Elementary') or np.any(np.array(self.reversible) =='yes'):
            raise NotImplementedError('Only elementary, irreversible reactions accepted')
        else:
            return self.reaction_rate(x)
