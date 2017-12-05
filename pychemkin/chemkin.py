from os import path
import numpy as np
from .SQLParser import SQLParser
from .InputParser import InputParser
from .ReactionCoeffs import ReactionCoeffs
from .BackwardCoeffs import BackwardCoeffs


class chemkin:

    '''
    The chemkin module computes the the reaction rates/ progress rates of all species
    given species concentrations and the systm temperature

    METHODS and ATTRIBUTES
    ========
    After initialization, user could call:
     - reaction_rate(x, T): method to calculate reaction rate given concentration x and temperature T.
     - species: A sequence (list, array, tuple) containing the names of the species
     - progress_rate(x): calculate progress rate given x (assumes temperature is already set)
     - _set_rc_params(T=..., R=..., A=...): internal method to set params of reaction coeffs

    EXAMPLES
    =========
    >>> chem = chemkin("tests/test_xml/rxns.xml")
    Finished reading xml input file
    >>> print(chem.species)
    ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    >>> chem.reaction_rate([1,1,1,1,1,1], 1000)
    array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
            -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    '''

    def __init__(self, file_name):
        '''
        INPUT
        =====
        file_name: string, required
                   xml file containing species and chemical reaction information
        '''
        self.file_name = file_name

        input_ = InputParser(file_name)

        self.rxndata = input_.reactions
        rc_list = [ReactionCoeffs(**params) for params in input_.rate_coeff_params]
        equationlist = []
        rxn_types = []
        reversible = []

        for i, reaction in enumerate(self.rxndata):
            equationlist.append(reaction['equation'])
            rxn_types.append(reaction['type'])
            reversible.append(reaction['reversible'].strip().lower())

        reversible = np.array(reversible)

        # check whether `reversible` is valid
        if not np.all(np.logical_or(reversible == 'yes', reversible == 'no')):
            raise ValueError('`reversible` should be either "yes" or "no".')

        reversible = reversible == 'yes'
        if np.any(reversible):
            species_reversible = np.sum(input_.nu_react[:, reversible] + input_.nu_prod[:, reversible], axis=1) > 0
        else:
            species_reversible = np.array([False for _ in range(len(input_.species))])

        if np.any(species_reversible) and np.any(reversible):
            the_shape = (np.sum(species_reversible),np.sum(reversible))
            backward_coeffs = BackwardCoeffs(input_.nu_react[:, reversible][species_reversible,:].reshape(the_shape), \
                                             input_.nu_prod[:, reversible][species_reversible,:].reshape(the_shape), \
                                             np.array(input_.species)[species_reversible])
        else:
            backward_coeffs = BackwardCoeffs(np.array([]), np.array([]), np.array([]))

        self.nu_react = input_.nu_react
        self.nu_prod = input_.nu_prod

        if self.nu_prod.shape != self.nu_react.shape or len(rc_list) != self.nu_prod.shape[1]:
            raise ValueError("Dimensions not consistent!")
        if np.any(self.nu_react<0.0) or np.any(self.nu_prod<0.0):
            raise ValueError("Negative stoichiometric coefficients are prohibited!")

        self.rc_list = rc_list
        self.bc = backward_coeffs
        self.species = input_.species
        self.equations = equationlist
        self.rxn_types = rxn_types

        self.reversible = reversible
        self.species_reversible = species_reversible
        if np.any(np.array(self.rxn_types)!='Elementary'):
            raise NotImplementedError('Only elementary reactions accepted')
        self.T = None

    def _set_rc_params(self,**kwargs):
        ''' add new or change parameters for all reaction coefficients
        '''
        if 'T' in kwargs:
            self.T = kwargs['T']
        for rc in self.rc_list:
            rc.set_params(**kwargs)

    def __repr__(self):
        class_name = type(self).__name__
        args = "'{}'".format(self.file_name)
        return class_name + "(" + args +")"

    def __len__(self):
        '''return number of reactions'''
        return len(self.rc_list)

    def summary(self,brief=True):
        '''
        Return information about reaction system in a more human-readable format than the input xml file

        INPUT
        =====
        brief: Boolean, optional
               If set to true, only return the equations, species, and reversibility
        '''
        eqn_str = "chemical equations:\n[\n" + "\n".join([str(eq_) for eq_ in self.equations]) + "\n]"
        species_str = "species: " + str(self.species)
        reversible_str = "reversible: " + str(self.reversible)
        if brief:
            return "\n".join([eqn_str,species_str,reversible_str])
        nu_react_str = "nu_react:\n" + str(self.nu_react)
        nu_prod_str = "nu_prod:\n" + str(self.nu_prod)
        rc_str = "reaction coefficients:\n[\n" + "\n".join([str(rc_) for rc_ in self.rc_list]) + "\n]"
        rxn_str = "reaction types: " + str(self.rxn_types)
        return "\n".join([eqn_str, species_str,nu_react_str,nu_prod_str,rc_str, rxn_str, reversible_str])

    def __str__(self):
        return self.summary(brief=False)

    def _progress_rate_default_T(self,x,kf,kb):
        '''Return progress rate for a system of M reactions involving N species, with temperature having been stored as a class attribute
        This is a helper function for the ODE solver that is not intended to be called directly. progress_rate is the user-facing function that can compute progress rates

        INPUTS
        =======
        x: array or list of floats, required
            A length-N vector specifying concentrations of the N species
        kf: array of floats, required
            Length-N array of forward reaction rate coefficients
        kb: array of floats, required
            Length-N array of backward reaction rate coefficients

        RETURNS
        ========
        Length-N array of progress rates of species
        '''
        if self.T == None:
            raise ValueError("T has not been set! Please use self._set_rc_params() to set T first!")

        pr = kf * np.product(x ** self.nu_react, axis=0)

        if np.any(self.reversible):
            pr[self.reversible] = pr[self.reversible] -  kb * np.product(x ** self.nu_prod[:, self.reversible], axis=0)
        return pr

    def _progress_rate_init(self,x,T):
        '''This is a helper function that reshapes the concentration vector and computes the reaction rate coefficients in order to allow progress rates to be calculated,d.

        INPUTS
        =======
        x: array or list, required
           A length-N vector specifying concentrations of the N species
        T: float, required
           Temperature in K

        RETURNS
        =======
        x: concentration vector reshaped appropriately for vectorized operations
        kf: Length-N array of forward reaction rate coefficients
        kb: Length-N array of backward reaction rate coefficients
        '''
        self._set_rc_params(T=T)
        x = np.array(x).reshape(-1,1)
        if len(x) != self.nu_prod.shape[0]:
            raise ValueError("ERROR: The concentration vector x must be of length N, where N is the \
            number of species")
        #check that concentrations are all non-negative:
        if np.any(x<0):
            raise ValueError("ERROR: All the species concentrations must be non-negative")
        #initialize progress rate vector
        kf = np.array([rc.k_forward() for rc in self.rc_list])
        kb = self.bc.backward_coeffs(kf[self.reversible], self.T)
        return x,kf,kb

    def progress_rate(self, x, T):
        '''Return progress rate for a system of M reactions involving N species

        INPUTS
        =======
        x: array or list, required
           A length-N vector specifying concentrations of the N species
        T: float, required
           Temperature in K

        RETURNS
        ========
        Length-N array of reaction rates of species
        '''
        x, kf, kb = self._progress_rate_init(x,T)
        return self._progress_rate_default_T(x,kf,kb)

    def reaction_rate(self, x, T):
        '''
        Return reaction rates for a system of M reactions involving N species

        INPUTS
        =======
        x: array or list, required
           A length-N vector specifying concentration of each specie
        T: float, optional
           Temperature in K

        RETURNS
        ========
        R: Array of reaction rates of species (length N)
        '''
        r = self.progress_rate(x, T)

        return np.sum(r * (self.nu_prod-self.nu_react), axis=1)
