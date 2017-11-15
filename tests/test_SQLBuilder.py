from pychemkin import *


def test_from_txt():
    '''
    Test building SQL database from TXT file
    '''
    sql = SQLBuilder('tests/test_data/thermo.txt', save_xml=True)
    print(sql.sql2pandas('tests/test_data/thermo.sqlite'))
    
def test_from_xml():
    '''
    Test building SQL database from XML file
    '''
    sql = SQLBuilder('tests/test_data/thermo.xml')