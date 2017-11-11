'''Thermodynamics and Thermochemistry for Chemical Kinetics
This module contains a BackwardCoeffs class with methods for 
computing the backward reaction rates for a set of 
reversible, elementary reactions.
'''

import numpy as np

class BackwardCoeffs:
    '''Methods for calculating the backward reaction rate coefficients.
    Cp_over_R: Returns specific heat of each specie given by 
               the NASA polynomials.
    H_over_RT:  Returns the enthalpy of each specie given by 
                the NASA polynomials.
    S_over_R: Returns the entropy of each specie given by 
              the NASA polynomials.
    backward_coeffs:  Returns the backward reaction rate 
                      coefficient for reach reaction.
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
        self.gamma = np.sum(self.rxnset.nuij, axis=0)

    def Cp_over_R(self, T):
        
        if self.coeffs is None:
            self.coeffs = self.sql.get_multi_coeffs(self.species, T)
        a = self.coeffs

        Cp_R = (a[:,0] + a[:,1] * T + a[:,2] * T**2.0 
                + a[:,3] * T**3.0 + a[:,4] * T**4.0)

        return Cp_R

    def H_over_RT(self, T):
        
        if self.coeffs is None:
            self.coeffs = self.sql.get_multi_coeffs(self.species, T)
        a = self.coeffs

        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
                + a[:,5] / T)

        return H_RT
               

    def S_over_R(self, T):

        if self.coeffs is None:
            self.coeffs = self.sql.get_multi_coeffs(self.species, T)
        a = self.coeffs

        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

        return S_R

    def backward_coeffs(self, kf, T):
        
        self.coeffs = self.sql.get_multi_coeffs(self.species, T)

        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(self.nu.T, self.H_over_RT(T))
        delta_S_over_R = np.dot(self.nu.T, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction 
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = self.p0 / self.R / T

        # Ke
        kb = fact**self.gamma * np.exp(delta_G_over_RT)

        return kf / kb