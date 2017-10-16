import xml.etree.ElementTree as ET
import numpy as np
from copy import deepcopy

class InputParser:
    """
    This class Input Parser takes the input of a xml file 
    and creates a dictionary. 
   
    
    INPUTS
    =======
  'equation' - an equation of the chemical reaction to solve
  'id'- the number of reactions in the file
  'products'- the output of the chemical equation 
  'rateCoeffParams'- the variables needed to calculate k, or k
  'reactants' - the input of the chemical equation 
  'reversible'- yes/no if a reversable equation
  'type'- the type of reaction (i.e. 'elementary')

    RETURNS
    ========
   ###### returns a dictionary with the following keys: 
  'equation'
  'id'
  'products'
  'rateCoeffParams'
  'reactants'
  'reversible': 
  'type'
    
    EXAMPLES
    =========
    >>> input_ = InputParser('rxns.xml')
    >>> print(input_.species)
    ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    """
    
    
    def __init__(self, file_name):
        """
        __init__ is used whenever an object of the class is constructed. 

        """ 
        self.file_name = file_name
        self.raw = ET.parse(self.file_name).getroot()
        self.species = self.raw.find('phase').find('speciesArray').text.strip().split()
        self.reactions = self.get_reactions(self.raw)
        self.nu_react, self.nu_prod = self.get_nu(self.reactions, self.species)
        self.rate_coeff_params = self.get_rate_coeff_params(self.reactions)
        
    def get_reactions(self, raw):
        """
        Identifies: 
        - elements in reaction and number of moles, including reactants and products
        - k, Depending on the input (as a given constant or calculated using Arrhenius/modified Arrhenius functions)
        
        """
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
        """getter for rate coefficients"""
        return [reaction['rateCoeffParams'] for reaction in reactions]
    
    def __repr__(self):
        '''Return a printable representation of the object.'''
        return 'InputParser(file_name=\'{}\')'.format(self.file_name)
    
    def __len__(self):
        '''Return the number of chemical reactions.'''
        return len(self.reactions)
    
