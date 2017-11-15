import xml.etree.ElementTree as ET
import re
from copy import deepcopy
from bs4 import BeautifulSoup
import numpy as np
import sqlite3
import pandas as pd

class SQLBuilder:
    '''
    Class for importing thermodynamic data from an XML file into an sqlite database
    '''
    def __init__(self, file_name, sql_name=None, save_xml=False):
        self.data = None
        if file_name[-4:] == '.txt':
            self.raw2sql(file_name, sql_name, save_xml)
        elif file_name[-4:] == '.xml':
            self.xml2sql(file_name, sql_name)

    def raw2data(self, file_name, save_xml=False, xml_name=None, i_start=5, lines=4):
        with open(file_name, 'r') as f:
            raw = f.read().strip()
        raw = raw.split('\n')
        species = [raw[i:i+lines] for i in range(i_start, len(raw), lines)][:-1]

        sci_number = re.compile('-?[0-9]+\.?[0-9]*E[+-]?[0-9]*')
        def find_sci_numbers(s):
            return [x for x in re.findall(sci_number, s)]

        def parse_specie(specie):
            name = specie[0].split()[0].strip()
            Ts = specie[0].split()[-4:-1]
            line_vals = [find_sci_numbers(s) for s in specie[1:]]
            coeff_high = line_vals[0] + line_vals[1][:2]
            coeff_low = line_vals[1][2:] + line_vals[2]
            return name, Ts, coeff_high, coeff_low

        data = dict()
        data['speciesArray'] = []
        data['speciesData'] = dict()

        for specie in species:
            name, Ts, coeff_high, coeff_low = parse_specie(specie)
            data['speciesArray'].append(name)
            data['speciesData'][name] = {'Ts':Ts, 'coeff_high':coeff_high, 'coeff_low':coeff_low}
        if save_xml:
            if xml_name is None:
                xml_name = file_name[:-4] + '.xml'
            self.rawdata2xml(data, xml_name)

        def rawdata2data(rawdata):
            data = deepcopy(rawdata)
            for name, _data in data['speciesData'].items():
                _data['Ts'] = [float(s) for s in _data['Ts']]
                _data['coeff_high'] = [float(s) for s in _data['coeff_high']]
                _data['coeff_low'] = [float(s) for s in _data['coeff_low']]
            return data

        self.data = rawdata2data(data)

        return self

    def rawdata2xml(self, data, file_name):
        root = ET.Element('ctml')

        root.append(ET.Comment('phase gri30'))

        phase = ET.SubElement(root, 'phase', id='gri30')
        ET.SubElement(phase, 'speciesArray', datasrc='#species_data').text = ' '.join(data['speciesArray'])

        root.append(ET.Comment('species definitions'))

        def add_speciesData(speciesData, data, p0="100000.0"):
            for name in data['speciesArray']:
                _data = data['speciesData'][name]
                speciesData.append(ET.Comment('species {}'.format(name)))
                specie = ET.SubElement(speciesData, 'species', name=name)
                thermo = ET.SubElement(specie, 'thermo')
                low = ET.SubElement(thermo, 'NASA', Tmax=_data['Ts'][-1], Tmin=_data['Ts'][0], p0=p0)
                ET.SubElement(low, 'floatArray', name='coeffs', size=str(len(_data['coeff_low']))).\
                text = ', '.join(_data['coeff_low'])
                high = ET.SubElement(thermo, 'NASA', Tmax=_data['Ts'][1], Tmin=_data['Ts'][-1], p0=p0)
                ET.SubElement(high, 'floatArray', name='coeffs', size=str(len(_data['coeff_high']))).\
                text = ', '.join(_data['coeff_high'])

        speciesData = ET.SubElement(root, 'speciesData', id='species_data')
        add_speciesData(speciesData, data)

        tree = ET.ElementTree(root)

        xmlstr = BeautifulSoup(ET.tostring(root), 'xml').prettify()
        with open(file_name, 'w') as f:
            f.write(xmlstr)

        return self


    def xml2data(self, file_name):
        raw = ET.parse(file_name).getroot()
        data = dict()
        data['speciesArray'] = raw.find('phase').find('speciesArray').text.strip().split()
        data['speciesData'] = dict()
        for specie in raw.find('speciesData').findall('species'):
            name = specie.attrib['name']
            data['speciesData'][name] = dict()
            data['speciesData'][name]['Ts'] = [0 for _ in range(3)]
            NASAs = specie.find('thermo').findall('NASA')
            i_high = np.argmax([float(NASA.attrib['Tmax']) for NASA in NASAs])
            data['speciesData'][name]['Ts'][0] = float(NASAs[1 - i_high].attrib['Tmin'])
            data['speciesData'][name]['Ts'][1] = float(NASAs[i_high].attrib['Tmax'])
            data['speciesData'][name]['Ts'][2] = float(NASAs[i_high].attrib['Tmin'])
            data['speciesData'][name]['coeff_high'] = [float(s) for s in NASAs[i_high].\
                                                       find('floatArray').text.strip().split(',')]
            data['speciesData'][name]['coeff_low'] = [float(s) for s in NASAs[1 - i_high].\
                                                       find('floatArray').text.strip().split(',')]
        self.data = data

        return self

    def data2sql(self, data, sql_name):
        db = sqlite3.connect(sql_name)
        cursor = db.cursor()
        cursor.execute('DROP TABLE IF EXISTS LOW')
        cursor.execute('DROP TABLE IF EXISTS HIGH')
        cursor.execute('''CREATE TABLE LOW (
                        SPECIES_NAME TEXT PRIMARY KEY NOT NULL,
                        TLOW FLOAT,
                        THIGH FLOAT,
                        COEFF_1 FLOAT,
                        COEFF_2 FLOAT,
                        COEFF_3 FLOAT,
                        COEFF_4 FLOAT,
                        COEFF_5 FLOAT,
                        COEFF_6 FLOAT,
                        COEFF_7 FLOAT)''')
        cursor.execute('''CREATE TABLE HIGH (
                        SPECIES_NAME TEXT PRIMARY KEY NOT NULL,
                        TLOW FLOAT,
                        THIGH FLOAT,
                        COEFF_1 FLOAT,
                        COEFF_2 FLOAT,
                        COEFF_3 FLOAT,
                        COEFF_4 FLOAT,
                        COEFF_5 FLOAT,
                        COEFF_6 FLOAT,
                        COEFF_7 FLOAT)''')
        for species in data['speciesArray']:
            _data = data['speciesData'][species]
            _low = [species, _data['Ts'][0], _data['Ts'][-1]] + _data['coeff_low']
            _high = [species, _data['Ts'][-1], _data['Ts'][1]] + _data['coeff_high']
            cursor.execute('''INSERT INTO LOW (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2,
            COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', _low)
            cursor.execute('''INSERT INTO HIGH (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2,
            COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', _high)
        db.commit()
        db.close()
        return self

    def raw2sql(self, raw_name, sql_name=None, save_xml=False):
        if sql_name is None:
            sql_name = raw_name[:-4] + '.sqlite'
        self.raw2data(raw_name, save_xml)
        self.data2sql(self.data, sql_name)
        return self

    def xml2sql(self, xml_name, sql_name=None):
        if sql_name is None:
            sql_name = xml_name[:-4] + '.sqlite'
        self.xml2data(xml_name)
        self.data2sql(self.data, sql_name)
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
