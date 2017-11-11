import numpy as np
import sqlite3
import pandas as pd

class SQLParser:
    def __init__(self, sql_name):
        self.data = None
        self.sql2data(sql_name)
    
    def sql2data(self, sql_name):
        db = sqlite3.connect(sql_name)
        cursor = db.cursor()
        data = dict()
        for row in cursor.execute('''SELECT * FROM LOW''').fetchall():
            data[row[0]] = {'low':{'Ts':row[1:3], 'coeffs':row[3:]}}
        for row in cursor.execute('''SELECT * FROM HIGH''').fetchall():
            if row[0] not in data:
                data[row[0]] = {'high':{'Ts':row[1:3], 'coeffs':row[3:]}}
            else:
                data[row[0]]['high'] = {'Ts':row[1:3], 'coeffs':row[3:]}
        db.close()
        self.data = data
        return self
    
    def sql2pandas(self, sql_name):
        db = sqlite3.connect(sql_name)
        cursor = db.cursor()
        cols = ['SPECIES_NAME', 'TLOW', 'THIGH', 'COEFF_1', 'COEFF_2', 'COEFF_3', \
                'COEFF_4', 'COEFF_5', 'COEFF_6', 'COEFF_7']
        queries = ['''SELECT * FROM LOW''', '''SELECT * FROM HIGH''']
        qs = [cursor.execute(query).fetchall() for query in queries]
        dfs = [pd.DataFrame.from_items([(col_name, [col[i] for col in q]) \
                                        for i, col_name in enumerate(cols)]) for q in qs]
        db.close()
        return dfs
    
    def get_coeffs(self, species, temp):
        '''Return coeffs of species at temp'''
        if not species in self.data:
            raise ValueError('Species not found in the provided NASA polynomials database.')
        data = self.data[species]
        if 'low' in data and data['low']['Ts'][0] <= temp and temp <= data['low']['Ts'][1]:
            return data['low']['coeffs']
        elif 'high' in data and data['high']['Ts'][0] <= temp and temp <= data['high']['Ts'][1]:
            return data['high']['coeffs']
        else:
            raise ValueError('Temperature not supported for the species in the \
            provided NASA polynomials database.')
    
    def get_species(self, temp):
        '''Return a list of supported species at temp'''
        species_list = []
        for species, data in self.data.items():
            if ('low' in data and data['low']['Ts'][0] <= temp and temp <= data['low']['Ts'][1]) \
            or ('high' in data and data['high']['Ts'][0] <= temp and temp <= data['high']['Ts'][1]):
                species_list.append(species)
        return species_list