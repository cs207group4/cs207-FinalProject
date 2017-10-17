from chemkin import *
def test_not_implement():
    """Test that input parser correctly handles reaction types and rate coefficients that are not implemented"""
    try:
        reactions = chemkin.from_xml('test_xml/rxns_bad_not_implement.xml')
    except NotImplementedError as e:
        assert type(e) == NotImplementedError
        print(e)
    try:
        reactions = chemkin.from_xml('test_xml/rxns_bad_not_implement_reverse.xml')
    except NotImplementedError as e:
        assert type(e) == NotImplementedError
        print(e)

def test_bad_reactant():
    """test that input parser correctly handles reactants not matching the species array"""
    try:
        reactions = chemkin.from_xml('test_xml/rxns_bad_reactants.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_correct():
    """test of standard usage"""
    input_ = InputParser('test_xml/rxns.xml')
    assert len(input_) == 3
    assert len(eval(repr(input_))) == 3
    print(repr(input_))

def test_singlerxn():
    """test that input parser correctly handles the single reaction case"""
    input_ = InputParser('test_xml/rxns_single.xml')
    assert len(input_) == 1
    assert np.array_equal(input_.nu_prod, np.array([[0],[0], [1], [1]]))
    assert np.array_equal(input_.nu_react, np.array([[1],[1], [0], [0]]))
    assert input_.rate_coeff_params[0]['k'] == 10000.0
    assert input_.rate_coeff_params[0]['type'] == 'Constant'
    print(repr(input_))

def test_extradata():
    """test that input parser correctly identifies too many reaction systems present"""
    try:
        input_ = InputParser('test_xml/rxns_doubleddata.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_missingspecies():
    """test that input parser correctly identifies that species array is missing"""
    try:
        input_ = InputParser('test_xml/rxns_missingspecies.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)

def test_missingcoeffdata():
    """test that input parser correctly identifies that coefficient data is missing"""
    try:
        input_ = InputParser('test_xml/rxns_missingcoeffdata.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)
