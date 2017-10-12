import xml.etree.ElementTree as ET
import numpy as np
from copy import deepcopy

class InputParser:
    
    def __init__(self, file_name):
        self.raw = ET.parse(file_name).getroot()
        self.species = self.raw.find('phase').find('speciesArray').text.split()
        self.reactions = self.get_reactions(self.raw)
        self.nu_react, self.nu_prod = self.get_nu(self.reactions, self.species)
        self.rate_coeff_params = self.get_rate_coeff_params(self.reactions)
        
    def get_reactions(self, raw):
        
        def parse_rate_coeff(reaction, reaction_dict):
            rc_ = reaction.find('rateCoeff')
            reaction_dict['rateCoeffParams'] = dict()
            if None != rc_.find('Constant'):
                reaction_dict['rateCoeffParams']['type'] = 'Constant'
                reaction_dict['rateCoeffParams']['k'] = float(rc_.find('Constant').find('k').text)
            elif None != rc_.find('Arrhenius'):
                reaction_dict['rateCoeffParams']['type'] = 'Arrhenius'
                reaction_dict['rateCoeffParams']['A'] = float(rc_.find('Arrhenius').find('A').text)
                reaction_dict['rateCoeffParams']['E'] = float(rc_.find('Arrhenius').find('E').text)
            elif None != rc_.find('modifiedArrhenius'):
                reaction_dict['rateCoeffParams']['type'] = 'modifiedArrhenius'
                reaction_dict['rateCoeffParams']['A'] = float(rc_.find('modifiedArrhenius').find('A').text)
                reaction_dict['rateCoeffParams']['b'] = float(rc_.find('modifiedArrhenius').find('b').text)
                reaction_dict['rateCoeffParams']['E'] = float(rc_.find('modifiedArrhenius').find('E').text)
            else:
                raise NotImplementedError('The type of reaction rate coefficient has not been implemented.')
            
        reactions = []
        for i, reaction in enumerate(raw.find('reactionData')):
            reactions.append(deepcopy(reaction.attrib))
            parse_rate_coeff(reaction, reactions[i])
            reactions[i]['equation'] = reaction.find('equation').text
            reactions[i]['reactants'] = {s.split(':')[0]:float(s.split(':')[1]) \
                                         for s in reaction.find('reactants').text.split()}
            reactions[i]['products'] = {s.split(':')[0]:float(s.split(':')[1]) \
                                         for s in reaction.find('products').text.split()}
        return reactions
    
    def get_nu(self, reactions, species):
        nu_react = np.zeros((len(species), len(reactions)))
        nu_prod = np.zeros((len(species), len(reactions)))
        for i, reaction in enumerate(reactions):
            if not (reaction['reversible'] == 'no' and reaction['type'] == 'Elementary'):
                raise NotImplementedError('The type of reaction has not been implemented.')
            for specie, stoi in reaction['reactants'].items():
                nu_react[species.index(specie), i] = stoi
            for specie, stoi in reaction['products'].items():
                nu_prod[species.index(specie), i] = stoi
        return nu_react, nu_prod
    
    def get_rate_coeff_params(self, reactions):
        return [reaction['rateCoeffParams'] for reaction in reactions]import numpy as np
class reaction_coeffs:
    ''' A class for reaction_coeffs
        
        Initialize by calling reaction_coeffs(type, {params})
        Set params by calling set_params({params})
        Get calculated k value by calling kval()
        
        EXAMPLES
        =========
        >>> rc = reaction_coeffs('Constant', k = 1e3)
        >>> rc.kval() 
        1000.0
        >>> rc = reaction_coeffs('Arrhenius', A = 1e7, E=1e3)
        >>> rc.set_params(T=1e2)
        >>> rc.kval() 
        3003549.0889639612
    '''

    def __init__(self, rtype, **kwargs):
        if rtype not in ["Constant","Arrhenius","modifiedArrhenius"]:
            raise NotImplementedError('Type Not Supported Yet!')
        self.__rtype = str(rtype)
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
        
        EXAMPLES
        =========
        >>> __const(1.0)
        1.0
        
        '''
        
        if k < 0:
            raise ValueError("k can not be negative!")
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
        
        EXAMPLES
        =========
        >>> __arr(1e7,1e3,1e2)
        3003549.0889639612
        
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
        
        EXAMPLES
        =========
        >>> __mod_arr(1e7,0.5,1e3,1e2)
        30035490.889639609
        
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
            return np.inf