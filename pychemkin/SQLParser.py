import os
import numpy as np
import sqlite3
import pandas as pd

class SQLParser:
    '''Class to retrieve relevant NASA thermodynamic data from sqlite database

    _sql2data: Query SQL database and store NASA polynomial coefficients in a dictionary
    _get_coeffs: Helper function to return NASA polynomial coefficients of a specie at temperature T
    get_multi_coeffs: Return NASA polynomial coefficients for an array of species at temperature T
    get_species: Return a list of species with valid NASA polynomial coefficients provided for temperature T

    EXAMPLES
    =========
    >>> parser = SQLParser()
    >>> parser.get_multi_coeffs(['O','O2'], 700)
    array([[  3.16826710e+00,  -3.27931884e-03,   6.64306396e-06,
             -6.12806624e-09,   2.11265971e-12,   2.91222592e+04,
              2.05193346e+00],
           [  3.78245636e+00,  -2.99673416e-03,   9.84730201e-06,
             -9.68129509e-09,   3.24372837e-12,  -1.06394356e+03,
              3.65767573e+00]])

    '''
    def __init__(self):
        '''
        INPUT
        =====
        sql_name: string, required
                  Name of sqlite database containing the NASA thermodynamic data
        '''
        self.data = None
        here = os.path.abspath(os.path.dirname(__file__))
        self.sql_name = os.path.join(here, 'data/thermo_all.sqlite')
        if not os.path.isfile(self.sql_name):
            raise ValueError('Database {} does not exist'.format(self.sql_name))
        self._sql2data()

    def _sql2data(self):
        '''Query SQL database and store NASA polynomial coefficients in a dictionary
        '''
        db = sqlite3.connect(self.sql_name)
        cursor = db.cursor()
        data = dict()
        #get low-temperature coefficients
        for row in cursor.execute('''SELECT * FROM LOW''').fetchall():
            data[row[0]] = {'low':{'Ts':row[1:3], 'coeffs':row[3:]}}
        #get high-temperature coefficients
        for row in cursor.execute('''SELECT * FROM HIGH''').fetchall():
            if row[0] not in data:
                data[row[0]] = {'high':{'Ts':row[1:3], 'coeffs':row[3:]}}
            else:
                data[row[0]]['high'] = {'Ts':row[1:3], 'coeffs':row[3:]}
        db.close()
        self.data = data


    def _get_coeffs(self, specie, T):
        '''Helper function to return NASA polynomial coefficients of a specie at temperature T

        INPUT
        =====
        specie: string, required
                specie for which NASA polynomial coefficients are desired
        T: float, required
           Temperature of reactor

        RETURNS
        =======
        length-7 tuple of NASA polynomial coefficients for a given specie at a given temperature

        '''
        if not specie in self.data:
            raise ValueError('Specie {} not found in the provided NASA polynomials database.'.format(specie))
        data = self.data[specie]
        if 'low' in data and data['low']['Ts'][0] <= T and T <= data['low']['Ts'][1]:
            return data['low']['coeffs']
        elif 'high' in data and data['high']['Ts'][0] < T and T <= data['high']['Ts'][1]:
            return data['high']['coeffs']
        else:
            raise ValueError\
            ('{} K out of valid temperature range for specie {} in the provided NASA polynomials database.'.format(T, specie))

    def get_multi_coeffs(self, species_array, T):
        '''Return NASA polynomial coefficients for an array of species at temperature T

        INPUT
        =====
        species_array: array of strings, required
                       Array of N species for which NASA polynomial coefficients are desired
        T: float, required
           Temperature of reactor

        RETURNS
        ======
        Nx7 array where each row corresponds to the NASA polynomial coefficients for the corresponding specie

        '''
        return np.array([self._get_coeffs(specie, T) for specie in species_array])
    #does this function get used anywhere?
    def get_species(self, T):
        '''Return a list of species with valid NASA polynomial coefficients provided for temperature T
        INPUT
        =====
        T: float, required
           Temperature of reactor

        RETURNS
        =======
        List of species with valid NASA polynomial coefficients at the requested temperature
        '''
        species_list = []
        for species, data in self.data.items():
            if ('low' in data and data['low']['Ts'][0] <= T and T <= data['low']['Ts'][1]) \
            or ('high' in data and data['high']['Ts'][0] < T and T <= data['high']['Ts'][1]):
                species_list.append(species)
        return species_list


    def sql2pandas(self):
        '''
        '''
        db = sqlite3.connect(self.sql_name)
        cursor = db.cursor()
        cols = ['SPECIES_NAME', 'TLOW', 'THIGH', 'COEFF_1', 'COEFF_2', 'COEFF_3', \
                'COEFF_4', 'COEFF_5', 'COEFF_6', 'COEFF_7']
        queries = ['''SELECT * FROM LOW''', '''SELECT * FROM HIGH''']
        qs = [cursor.execute(query).fetchall() for query in queries]
        dfs = [pd.DataFrame.from_items([(col_name, [col[i] for col in q]) \
                                        for i, col_name in enumerate(cols)]) for q in qs]
        db.close()
        return dfs
