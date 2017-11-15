import numpy as np
from .InputParser import InputParser
from .SQLParser import SQLParser

class BackwardCoeffs:
    '''Methods for calculating the backward reaction rate coefficients for a set of reversible, elementary reactions.
    Cp_over_R: Return specific heat of each specie given by
               the NASA polynomials.
    H_over_RT:  Return the enthalpy of each specie given by
                the NASA polynomials.
    S_over_R: Return the entropy of each specie given by
              the NASA polynomials.
    backward_coeffs:  Return the backward reaction rate
                      coefficient for each reaction.

    EXAMPLE
    =======
    >>> ip = InputParser('tests/test_xml/rxns_reversible.xml')
    Finished reading xml input file
    >>> sql = SQLParser('pychemkin/data/thermo.sqlite')
    >>> bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    '''

    def __init__(self, nu_react, nu_prod, species, sql):
        '''
        INPUT
        =====
        nu_prod: Array of integers, required
                 N X M array of stoichiometric coefficients for products (N species, M reactions)
        nu_react: Array of integers, required
                 N x M array of stoichiometric coefficients for reactants
        species: Array of strings, required
                 List or array of length N providing the species
        sql: Instance of SQLParser, required

        '''
        self.nu = nu_prod - nu_react
        self.species = species
        self.sql = sql
        self.coeffs = None

        self.p0 = 1.0e+05 # Pa
        self.R = 8.3144598 # J / mol / K
        self.gamma = np.sum(self.nu, axis=0)

    def Cp_over_R(self, T):
        '''
        Return specific heat of each specie given by the NASA polynomials.

        INPUT
        =====
        T: float, required
           Temperature of reactor

        RETURNS
        =======
        Float or array of floats corresponding to the heat capacity at constant pressure divided by the ideal constant
        '''

        if self.coeffs is None:
            self.coeffs = self.sql.get_multi_coeffs(self.species, T)
        a = self.coeffs

        Cp_R = (a[:,0] + a[:,1] * T + a[:,2] * T**2.0
                + a[:,3] * T**3.0 + a[:,4] * T**4.0)

        return Cp_R

    def H_over_RT(self, T):
        '''
        Return enthalpy of each specie given by the NASA polynomials.

        INPUT
        =====
        T: float, required
           Temperature of reactor

        RETURNS
        =======
        Float or array of floats corresponding to enthalpy divded by RT

        '''

        if self.coeffs is None:
            self.coeffs = self.sql.get_multi_coeffs(self.species, T)
        a = self.coeffs

        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0
                + a[:,5] / T)

        return H_RT


    def S_over_R(self, T):
        '''
        Return entropy of each specie given by the NASA polynomials.

        INPUT
        =====
        T: float, required
           Temperature of reactor

        RETURNS
        =======
        Float or array of floats corresponding to entropy divided by the ideal gas constant

        '''

        if self.coeffs is None:
            self.coeffs = self.sql.get_multi_coeffs(self.species, T)
        a = self.coeffs

        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

        return S_R

    def backward_coeffs(self, kf, T):
        '''
        Calculate backward reaction rate coefficients for a system of M reactions
        INPUT
        =====
        kf: float or float array, required
            forward reaction rate coefficients
        T: float, required
           Temperature of reactor
        '''

        if len(kf) == 0:
            return np.array([])

        self.coeffs = self.sql.get_multi_coeffs(self.species, T)

        if len(kf)!=len(self.gamma):
            raise ValueError('Array of forward reaction rate coefficients needs to match with ')


        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(self.nu.T, self.H_over_RT(T))
        delta_S_over_R = np.dot(self.nu.T, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = self.p0 / self.R / T

        # Ke
        ke = fact**self.gamma * np.exp(delta_G_over_RT)

        return kf / ke
