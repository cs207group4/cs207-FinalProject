import xml.etree.ElementTree as ETs
import numpy as np
from collections import OrderedDict
from .SQLParser import SQLParser
from .InputParser import InputParser
from .ReactionCoeffs import ReactionCoeffs
from .BackwardCoeffs import BackwardCoeffs
from os import path

class chemkin:

    '''
    This class Chemkin computes the the reaction rates/ progress rates of the species
    at each temperature of interest given species concentrations
    INPUTS
    =======
    Using an XML file, the class calls InputParser and ReactionCoeffs to calculate the reaction rates
    RETURNS
    ========
    After initialization, user could call:
     - reaction_rate(x, T): method to calculate reaction rate given concentration x and temperature T.
     - species: A sequence (list, array, tuple) containing sthe names of the species
     - progress_rate(x): calculate progress rate given x (assumes temperature is already set)
     - _set_rc_params(T=..., R=..., A=...): internal method to set params of reaction coeffs
    EXAMPLES
    =========
    >>> chem = chemkin("tests/test_xml/rxns.xml", sql_name = "pychemkin/data/thermo.sqlite")
    Finished reading xml input file
    >>> print(chem.species)
    ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    >>> chem.reaction_rate([1,1,1,1,1,1], 1000)
    array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
            -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    '''

    def __init__(self, file_name, sql_name = None):
        '''
        INPUT
        =====
        file_name: string, required
                   xml file containing species and chemical reaction information
        sql_name: string, optional
                Name of sqlite database holding the NASA thermodynamic data
        '''
        self.file_name = file_name
        
        if sql_name==None:
            here = path.abspath(path.dirname(__file__))
            sql_name = path.join(here, 'data/thermo.sqlite')

        self.sql_name = sql_name
        sql = SQLParser(sql_name)
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
        reversible = np.array(reversible) == 'yes'

        backward_coeffs = BackwardCoeffs(input_.nu_react[:, reversible], input_.nu_prod[:, reversible], input_.species, sql)

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
        if np.any(np.array(self.rxn_types)!='Elementary'):
            raise NotImplementedError('Only elementary reactions accepted')
        self.T = None

    def _set_rc_params(self,**kwargs):
        ''' add new or change old parameters for all reaction coeffs.
        '''
        if 'T' in kwargs:
            self.T = kwargs['T']
        for rc in self.rc_list:
            rc.set_params(**kwargs)

    def __repr__(self):
        class_name = type(self).__name__
        args = "'{}', sql_name = '{}'".format(self.file_name, self.sql_name)
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

    def progress_rate(self, x, T):
        '''Return progress rate for a system of M reactions involving N species

        INPUTS
        =======
        x: array or list, required
           A length-N vector specifying concentrations of the N species
        T: float, optional
           Temperature in K
        RETURNS
        ========
        Length-N array of reaction rates of species
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
        pr = kf* np.product(x ** self.nu_react, axis=0)

        if np.any(self.reversible):
            pr[self.reversible] = pr[self.reversible] - self.bc.backward_coeffs(kf[self.reversible], self.T) \
            * np.product(x ** self.nu_prod[:, self.reversible], axis=0)
        return pr

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
        self._set_rc_params(T=T)
        r = self.progress_rate(x, T)

        # Return an array...
        return np.sum(r * (self.nu_prod-self.nu_react), axis=1)
