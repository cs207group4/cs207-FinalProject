import xml.etree.ElementTree as ET
from copy import deepcopy
import numpy as np


class InputParser:
    '''
    The Input Parser class processes xml input used for subsequent chemical kinetics calculations

    Attributes
    ========
    'equation' - an equation of the chemical reaction to solve
    'id'- the number of reactions in the file
    'products'- the output of the chemical equation
    'rateCoeffParams'- the variables needed to calculate the reaction rate coefficient
    'reactants' - chemical reactants
    'reversible'- yes/no if a reversible reaction
    'type'- the type of reaction (i.e. 'elementary')
    EXAMPLES
    =========
    >>> input_ = InputParser('tests/test_xml/rxns.xml')
    Finished reading xml input file
    >>> print(input_.species)
    ['H', 'O', 'OH', 'H2', 'H2O', 'O2']
    '''
    def __init__(self, file_name):
        '''
        INPUT
        =====
        file_name: string, required
                   name of xml input file
        '''
        self.file_name = file_name
        self.raw = ET.parse(self.file_name).getroot()
        self.species = self.get_species()
        self.reactions = self.get_reactions()
        self.nu_react, self.nu_prod = self.get_nu()
        self.rate_coeff_params = self.get_rate_coeff_params()
        print("Finished reading xml input file")

    def get_species(self):
        '''
        Returns species array from xml file
        '''
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
        '''
        Returns list of dictionaries of information about each reaction in the xml file
        '''
        def parse_rate_coeff(reaction, reaction_dict):
            rc_ = reaction.findall('rateCoeff')
            if len(rc_) > 1:
                raise ValueError('ERROR: there are multiple `rateCoeff`s under a reaction')
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
                    reactions[i]['reactants'] = {s.split(':')[0]:float(s.split(':')[1]) \
                                         for s in reaction.find('reactants').text.split()}
                    reactions[i]['products'] = {s.split(':')[0]:float(s.split(':')[1]) \
                                         for s in reaction.find('products').text.split()}

                except AttributeError:
                    print("Check that equation, reactants, and products information are present for every reaction")
                    raise

                except ValueError:
                    print("Stoichiometric coefficients must be float numbers")
                    raise
        return reactions

    def get_nu(self):
        '''
        Return tuple of arrays corresponding to stoichiometric coefficients for reactants and for products
        '''

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
