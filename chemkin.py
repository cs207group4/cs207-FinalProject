import xml.etree.ElementTree as ET
import numpy as np
from copy import deepcopy
from collections import OrderedDict

class InputParser:

    """
    This class Input Parser takes the input of a xml file
    and creates useful variables can be called by other functions or users


    INPUTS
    =======
    file_name: the xml filename.

    RETURNS
    ========
    ###### After init, following class variables can be used:
    'equation' - an equation of the chemical reaction to solve
    'id'- the number of reactions in the file
    'products'- the output of the chemical equation
    'rateCoeffParams'- the variables needed to calculate k, or k
    'reactants' - the input of the chemical equation
    'reversible'- yes/no if a reversable equation
    'type'- the type of reaction (i.e. 'elementary')

    EXAMPLES
    =========
    >>> input_ = InputParser('rxns.xml')
    Finished reading xml input file
    >>> print(input_.species)
    ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    """


    def __init__(self, file_name):
        """
        INPUT
        =====
        file_name: string, required
                   name of xml input file

        """
        self.file_name = file_name
        self.raw = ET.parse(self.file_name).getroot()
        self.species = self.get_species()
        self.reactions = self.get_reactions()
        self.nu_react, self.nu_prod = self.get_nu()
        self.rate_coeff_params = self.get_rate_coeff_params()
        print("Finished reading xml input file")

    def get_species(self):
        """
        Returns species array from xml file
        """
        try:
            species = self.raw.find('phase').find('speciesArray').text.strip().split()
        except AttributeError:
            print('ERROR: Check that valid species array is provided in xml input file')
            raise
        if len(species)==0:
            raise ValueError('ERROR: Species array needs to be provided')
        else:
            return species


    def get_reactions(self):
        """
        Returns list of dictionaries of information about each reaction in the xml file
        """
        def parse_rate_coeff(reaction, reaction_dict):
            rc_ = reaction.findall('rateCoeff')
            if len(rc_) > 1:
                raise NotImplementedError('ERROR: Only irreversible reactions are currently supported')
            elif len(rc_)==0:
                raise ValueError('Rate coefficient data appear to be missing')
            rc_ = rc_[0]
            reaction_dict['rateCoeffParams'] = dict()
            if rc_.find('Constant') is not None:
                reaction_dict['rateCoeffParams']['type'] = 'Constant'
                try:
                    k = float(rc_.find('Constant').find('k').text)

                except AttributeError:
                    print('ERROR: Value of k must be provided for constant rate coefficient')
                    raise
                reaction_dict['rateCoeffParams']['k'] = k

            elif rc_.find('Arrhenius') is not None:
                reaction_dict['rateCoeffParams']['type'] = 'Arrhenius'
                try:
                    A = float(rc_.find('Arrhenius').find('A').text)
                    E = float(rc_.find('Arrhenius').find('E').text)
                except AttributeError:
                    print('ERROR: Values of A and E must be provided for Arrhenius rate coefficient')
                    raise
                reaction_dict['rateCoeffParams']['A'] = A
                reaction_dict['rateCoeffParams']['E'] = E

            elif rc_.find('modifiedArrhenius') is not None:
                reaction_dict['rateCoeffParams']['type'] = 'modifiedArrhenius'
                try:
                    A= float(rc_.find('modifiedArrhenius').find('A').text)
                    b = float(rc_.find('modifiedArrhenius').find('b').text)
                    E = float(rc_.find('modifiedArrhenius').find('E').text)
                except AttributeError:
                    print('ERROR: Values of A, b, and E must be provided for modified Arrhenius rate coefficient')
                    raise
                reaction_dict['rateCoeffParams']['A'] = A
                reaction_dict['rateCoeffParams']['b'] = b
                reaction_dict['rateCoeffParams']['E'] = E
            else:
                raise NotImplementedError('This type of reaction rate coefficient has not been implemented. Current supported types are constant, Arrhenius, and modified Arrhenius.')

        reactions = []
        #check that only one system of reactions is present in the file
        rxndata = self.raw.findall('reactionData')
        if len(rxndata) > 1:
            raise ValueError('ERROR: Only one reaction system allowed in input file')
        elif len(rxndata)==0 or len(rxndata[0])==0:
            raise ValueError('Error: No reactions appear to be present in input file')

        else:
            for i, reaction in enumerate(rxndata[0]):
                reactions.append(deepcopy(reaction.attrib))
                parse_rate_coeff(reaction, reactions[i])
                try:
                    reactions[i]['equation'] = reaction.find('equation').text
                    reactions[i]['reactants'] = {s.split(':')[0]:int(s.split(':')[1]) \
                                         for s in reaction.find('reactants').text.split()}
                    reactions[i]['products'] = {s.split(':')[0]:int(s.split(':')[1]) \
                                         for s in reaction.find('products').text.split()}

                except AttributeError:
                    print("Check that equation, reactants, and products information are present for every reaction")
                    raise

                except ValueError:
                    print("Stoichiometric coefficients must be integers")
                    raise
        return reactions

    def get_nu(self):
        """
        Return tuple of arrays corresponding to stoichiometric coefficients for reactants and for products
        """

        nu_react = np.zeros((len(self.species), len(self.reactions)), dtype = int)
        nu_prod = np.zeros((len(self.species), len(self.reactions)), dtype = int)

        for i, reaction in enumerate(self.reactions):
            for specie, stoi in reaction['reactants'].items():
                nu_react[self.species.index(specie), i] = stoi
            for specie, stoi in reaction['products'].items():
                nu_prod[self.species.index(specie), i] = stoi

        return nu_react, nu_prod

    def get_rate_coeff_params(self):
        """getter for rate coefficients"""
        return [reaction['rateCoeffParams'] for reaction in self.reactions]

    def __repr__(self):
        """Return a printable representation of the object."""
        return 'InputParser(file_name=\'{}\')'.format(self.file_name)

    def __len__(self):
        """Return the number of chemical reactions."""
        return len(self.reactions)

import numpy as np

class ReactionCoeffs:
    """
    A class for computating reaction coefficients

    Initialize by calling reaction_coeffs(type, {params})
    Set params by calling set_params({params})
    Calculate k value by calling kval()

    EXAMPLES
    =========
    >>> rc = ReactionCoeffs('Constant', k = 1e3)
    >>> rc.kval()
    1000.0
    >>> rc = ReactionCoeffs('Arrhenius', A = 1e7, E=1e3)
    >>> rc.set_params(T=1e2)
    >>> rc.kval()
    3003549.0889639617
    """

    def __init__(self, type, **kwargs):
        """
        INPUT
        =====
        type: string, required
              Type of the reaction rate coefficient of interest (e.g. "Constant", "Arrhenius", etc)
        """
        self.__rtype = str(type)
        self.__params = kwargs
        if 'R' not in self.__params:
            self.__params['R'] = 8.314

    def set_params(self, **kwargs):
        for param in kwargs:
            self.__params[param] = kwargs[param]

    def __str__(self):
        return self.__rtype + " Reaction Coeffs: " + str(self.__params)

    def __repr__(self):
        class_name = type(self).__name__
        params_str = ""
        for param in self.__params:
            params_str += " ,"+param+" = "+str(self.__params[param])
        return str(class_name)+'("'+self.__rtype+'"'+params_str+')'

    def __eq__(self,other):
        """
        Check if two reaction rate coefficients are the same. They are same if the k values are equal.
        """
        return self.kval() == other.kval()

    def __check_param_in(self, param_list):
        """
        Returns dictionary of parameter values

        INPUT
        =====
        param_list: iterable, required
                    list of parameter namesof interest

        RETURNS
        =======
        Dictionary of parameters in provided list (returns None if mismatch present between parameter list and instance attributes)
        """
        param_dict = dict()
        for param in param_list:
            if param in self.__params:
                param_dict[param] = self.__params[param]
            else:
                return None
        return param_dict

    def kval(self):
        """
        Computes reaction coefficient

        NOTES
        ==========
        call corresponding private method to calculate k_values.
        Easy to extend... Just write a new private method for new type of coeffs.

        RETURNS
        ==========
        k: (real number) Reaction coefficient
        """
        if self.__rtype == "Constant":
            params = self.__check_param_in(['k'])
            if params != None:
                return self.__const(**params)
            else:
                raise ValueError('Insufficient Parameters')
        elif self.__rtype == "Arrhenius":
            params = self.__check_param_in(['A','E','T','R'])
            if params != None:
                return self.__arr(**params)
            else:
                raise ValueError('Insufficient Parameters')
        elif self.__rtype == "modifiedArrhenius":
            params = self.__check_param_in(['A','b','E','T','R'])
            if params != None:
                return self.__mod_arr(**params)
            else:
                raise ValueError('Insufficient Parameters')
        else:
            raise NotImplementedError("Type not supported yet!")

    def __const(self,k):
        """
        Return constant coefficient

        INPUTS
        =======
        k: float, required
           constant reaction rate coefficient

        RETURNS
        ========
        k: float
           constant reaction rate coefficient

        """

        if k < 0:
            raise ValueError("k cannot be negative!")
        return k

    def __arr(self,A,E,T,R):
        """
        Return Arrhenius reaction rate coefficient

        INPUTS
        =======
        A: float, required
           Arrhenius prefactor  A. A  is strictly positive
        E: float, required
           Activation energy
        T: float, required
           Temperature. T must be positive (assuming a Kelvin scale)
        R: float, required
           Ideal gas constant


        RETURNS
        ========
        k: float
           Arrhenius reaction rate coefficient

        NOTES
        ========
        R = 8.314 is the default ideal gas constant.
        It should be positive (except to convert units)

        """

        if A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

        return A * np.exp(-E / R / T)

    def __mod_arr(self,A,b,E,T,R):
        """
        Return modified Arrhenius reaction rate coefficient

        INPUTS
        =======
        A: float, required
           Arrhenius prefactor  A. A  is strictly positive
        b: float, required
           Modified Arrhenius parameter b. b must be real
        E: float, required
           Activation energy
        T: float, required
           Temperature. T must be positive (assuming a Kelvin scale)
        R: float, required
           The ideal gas constant


        RETURNS
        ========
        k: Modified Arrhenius reaction rate coefficient

        NOTES
        ========
        R=8.314 is the default ideal gas constant
        """
        if A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

        return A * (T**b) * np.exp(-E / R / T)


import numpy as np
from InputParser import InputParser
from SQLParser import SQLParser
from ReactionCoeffs import ReactionCoeffs
from BackwardCoeffs import BackwardCoeffs


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
     - reaction_rate(x): method to calculate reaction rate given concentration x (assumes temperature is 
       already set)
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

    def __init__(self, nu_react, nu_prod, reaction_coeffs, backward_coeffs, rxn_types, \
                 reversible, species, equations=None):
        """
        INPUT
        =====
        nu_prod: Array of integers, required
                 N X M array of stoichiometric coefficients for products (N species, M reactions)
        nu_react: Array of integers, required
                 N x M array of stoichiometric coefficients for reactants
        reaction_coeffs: Array, required
                 Length M list of reaction rate coefficients in the form of representations of ReactionCoeff 
                 objects
        backward_coeffs: Instance of BackwardCoeffs, required
        rxn_types: Array of strings, required
                 List or array of length M of reaction types
        reversible: Boolean array, required
                 List or array of length M denoting whether the reaction is reversible
        species: Array of strings, required
                 List or array of length N providing the species
        equations: Array of strings, optional
                 List or array of length M providing the reaction equations
        """
        self.nu_react = np.array(nu_react)
        self.nu_prod = np.array(nu_prod)
        if self.nu_prod.shape != self.nu_react.shape or len(reaction_coeffs) != self.nu_prod.shape[1]:
            raise ValueError("Dimensions not consistent!")
        if np.any(self.nu_react<0.0) or np.any(self.nu_prod<0.0):
            raise ValueError("Negative stoichiometric coefficients are prohibited!")
        self.rc_list = reaction_coeffs
        self.bc = backward_coeffs
        self.species = species
        self.equations = equations
        self.rxn_types = rxn_types
        if not np.array_equal(reversible, reversible.astype(bool)):
            raise TypeError('`reversible` should be a boolean array.')
        self.reversible = reversible
        if np.any(np.array(self.rxn_types)!='Elementary'):
            raise NotImplementedError('Only elementary reactions accepted')
        self.T = None

    @classmethod
    def from_xml(cls, file_name, sql_name='thermo30.sqlite'):
        """
        Alternate constructor starting with an xml input file

        INPUT
        =====
        file_name: string, required
                Name of input xml file
        sql_name: string, optional
                Name of input sqlite file
                 
        RETURNS
        =======
        initialized chemkin object
        """
        sql = SQLParser(sql_name)
        input_ = InputParser(file_name)
        rc_list = [ReactionCoeffs(**params) for params in input_.rate_coeff_params]
        rxndata = input_.reactions
        equationlist = []
        rxn_types = []
        reversible = []
        
        for i, reaction in enumerate(rxndata):
            equationlist.append(reaction['equation'])
            rxn_types.append(reaction['type'])
            reversible.append(reaction['reversible'].strip().lower())
        reversible = np.array(reversible) == 'yes'
        
        bc = BackwardCoeffs(input_.nu_react[:, reversible], input_.nu_prod[:, reversible], input_.species, sql)
        
        return cls(input_.nu_react, input_.nu_prod, rc_list, bc, rxn_types, reversible, species=input_.species, \
                   equations=equationlist)

    @classmethod
    def init_const_rc(cls, nu_react, nu_prod, rcs, rxn_types, reversible, species, equations=None, \
                     sql_name='thermo30.sqlite'):
        ''' construct an object with all constant or precalculated coeffs.
        '''
        sql= SQLParser(sql_name)
        rc_list = [ReactionCoeffs(type="Constant", k=rc_) for rc_ in rcs]
        reversible = reversible.astype(bool)
        bc = BackwardCoeffs(np.array(nu_react)[:, reversible], np.array(nu_prod)[:, reversible], species, sql)
        return cls(nu_react, nu_prod, rc_list, bc, rxn_types, reversible, species, equations)

    def set_rc_params(self,**kwargs):
        ''' add new or change old parameters for all reaction coeffs.
        '''
        if 'T' in kwargs:
            self.T = kwargs['T']
        for rc in self.rc_list:
            rc.set_params(**kwargs)

    def __repr__(self):
        class_name = type(self).__name__
        args = repr(self.nu_react.tolist()) + ',' + repr(self.nu_prod.tolist()) + ',' + repr(self.rc_list) \
        + ',' + repr(self.bc) + ',' + repr(self.rxn_types)+ ',' + repr(self.reversible)+',' + repr(self.species)\
        + ',' + repr(self.equations)
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

    def progress_rate(self, x):
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

        x = np.array(x).reshape(-1,1)
        if len(x) != self.nu_prod.shape[0]:
            raise ValueError("ERROR: The concentration vector x must be of length N, where N is the \
            number of species")
        #check that concentrations are all non-negative:
        if np.any(x<0):
            raise ValueError("ERROR: All the species concentrations must be non-negative")

        pr = np.array([rc.kval() for rc in self.rc_list]).astype(float) * np.product(x ** self.nu_react, axis=0)
        
        if np.sum(self.reversible) > 0:
            if self.T is None:
                raise ValueError('T not set.')
            pr[self.reversible] = pr[self.reversible] - self.bc.backward_coeffs(pr[self.reversible], self.T) \
            * np.product(x ** self.nu_prod[:, self.reversible], axis=0)
        return pr

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
        self.set_rc_params(T=T)
        return self.reaction_rate(x)