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
        file_name: string
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
        Retrieves species array from xml file
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
                if k<=0:
                    raise ValueError('k must be non-negative')
                else:
                    reaction_dict['rateCoeffParams']['k'] = k

            elif rc_.find('Arrhenius') is not None:
                reaction_dict['rateCoeffParams']['type'] = 'Arrhenius'
                try:
                    A = float(rc_.find('Arrhenius').find('A').text)
                    E = float(rc_.find('Arrhenius').find('E').text)
                except AttributeError:
                    print('ERROR: Values of A and E must be provided for Arrhenius rate coefficient')
                    raise
                if A <=0 or E<0:
                    raise ValueError('A must always be positive and E must be non-negative')
                else:
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
                if A <=0 or E<0 or np.iscomplex(b):
                    raise ValueError('ERROR: A must always be positive, E must be non-negative, and b must be real')
                else:
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
            if not (reaction['reversible'] == 'no' and reaction['type'] == 'Elementary'):
                raise NotImplementedError('Reactions that are reversible or non-elementary are not yet handled by this program.')
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


class ReactionCoeffs:
    ''' A class for reaction_coeffs

        Initialize by calling reaction_coeffs(type, {params})
        Set params by calling set_params({params})
        Get calculated k value by calling kval()

        EXAMPLES
        =========
        >>> rc = ReactionCoeffs('Constant', k = 1e3)
        >>> rc.kval()
        1000.0
        >>> rc = ReactionCoeffs('Arrhenius', A = 1e7, E=1e3)
        >>> rc.set_params(T=1e2)
        >>> rc.kval()
        3003549.0889639612
    '''

    def __init__(self, type, **kwargs):
        if type not in ["Constant","Arrhenius","modifiedArrhenius"]:
            raise NotImplementedError('Rate coefficient type not supported yet! Must be constant, Arrhenius, or modified Arrhenius.')
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
        ''' check if two coeffs are the same. They are same if k is the same...
        '''
        return self.kval() == other.kval()

    def __check_param_in(self, param_list):
        '''
        return None if not all params avaliable
        else return a dict containing all avaliable params
        '''
        param_dict = dict()
        for param in param_list:
            if param in self.__params:
                param_dict[param] = self.__params[param]
            else:
                return None
        return param_dict

    def kval(self):
        ''' A wrapper function

        NOTES
        ==========
        call corresponding private method to calculate k_values.
        Easy to extend... Just write a new private method for new type of coeffs.

        RETURNS
        ==========
        k: (real number) reaction coeffs.

        '''
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
        ''' return constant coefficient

        INPUTS
        =======
        k: constant coefficient

        RETURNS
        ========
        k constant coefficient


        '''

        if k < 0:
            raise ValueError("k cannot be negative!")
        return k

    def __arr(self,A,E,T,R):
        ''' return Arrhenius reaction rate coefficient

        INPUTS
        =======
        A: Arrhenius prefactor  A. A  is strictly positive
        E: Activation energy
        T: Temperature. T must be positive (assuming a Kelvin scale)
        R: Ideal gas constant


        RETURNS
        ========
        k: Arrhenius reaction rate coefficient

        NOTES
        ========
        R=8.314 is the ideal gas constant. It should never be changed (except to convert units)

        '''
        # R should never be changed
        if A <= 0:
            raise ValueError("A must be positive!")
        if T <= 0:
            raise ValueError("T must be positive!")

        # On overflow, return np.inf
        try:
            return A*np.exp(-E/(R*T))
        except OverflowError:
            print("WARNING: Overflow error")
            return np.inf

    def __mod_arr(self,A,b,E,T,R):
        ''' return modified Arrhenius reaction rate coefficient

        INPUTS
        =======
        A: Arrhenius prefactor  A. A  is strictly positive
        b: Modified Arrhenius parameter b. b must be real
        E: Activation energy
        T: Temperature. T must be positive (assuming a Kelvin scale)
        R: The ideal gas constant


        RETURNS
        ========
        k: Modified Arrhenius reaction rate coefficient

        NOTES
        ========
        R=8.314 is the ideal gas constant. It should never be changed (except to convert units)


        '''
        # R should never be changed
        if A <= 0:
            raise ValueError("A must be positive!")
        if T <= 0:
            raise ValueError("T must be positive!")
        if not isinstance(b, (float, int)):
            raise ValueError("b must be real!")

        # On overflow, return np.inf
        try:
            return A*(T**b)*np.exp(-E/(R*T))
        except OverflowError:
            print("WARNING: Overflow error")
            return np.inf


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
    def init_const_rc(cls,nu_react,nu_prod,rcs,species=None,equations=None):
        ''' construct an object with all constant or precalculated coeffs.
        '''
        rc_list = [ReactionCoeffs(type="Constant", k=rc_) for rc_ in rcs]
        return cls(nu_react,nu_prod,rc_list,species,equations)

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
        if np.any(x<0):
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

        r = self.progress_rate(x)

        # Return an array...
        return np.sum(r * (self.nu_prod-self.nu_react), axis=1)

    def reaction_rate_T(self, x, T):
        '''
        A function to easily calculate reaction rate based on x and T.
        '''

        if T < 0:
            raise ValueError("ERROR: Temperature cannot be negative")
        self.set_rc_params(T=T)
        return self.reaction_rate(x)