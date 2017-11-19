from pychemkin import SQLParser

def test_databasexistence():
    '''
    Test whether program correctly catches that database doesn't exist (or has incorrect path)
    '''
    try:
        parser = SQLParser()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_missingspecie():
    '''
    Test whether SQL parser correctly handles a requested specie missing from the database
    '''
    try:
        parser = SQLParser()
        parser._get_coeffs('CH3CN', 700)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_outofTrange():
    '''
    Test whether SQL parser correctly handles a temperature request that's out of the valid range
    '''
    try:
        parser = SQLParser()
        parser._get_coeffs('O2',4000)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_get_species():
    '''
    Test get_species function
    '''
    parser = SQLParser()
    print(parser.get_species(1000))


def test_sql2pandas():
    '''
    Test sql2pandas function
    '''
    parser = SQLParser()
    print(parser.sql2pandas())
