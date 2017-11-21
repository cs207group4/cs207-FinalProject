#change this once we rearrange directory structure
from pychemkin import chemkin
import numpy as np

def test_not_implement():
    """Test that input parser correctly handles reaction types and rate coefficients that are not implemented"""
    try:
        reactions = chemkin('tests/test_xml/rxns_bad_not_implement.xml')
    except NotImplementedError as e:
        assert type(e) == NotImplementedError
        print(e)

def test_bad_reactant():
    """test that input parser correctly handles reactants not matching the species array"""
    try:
        reactions = chemkin('tests/test_xml/rxns_bad_reactants.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_from_xml():
    """
    testing normal reaction rate calculations from xml input
    """
    reactions = chemkin('tests/test_xml/rxns.xml')
    expect = np.array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
        -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    assert all(reactions.reaction_rate([1, 1, 1, 1, 1, 1], 1000).astype(int) == expect.astype(int))
    return reactions


def test_singlerxnxml():

    """testing that single reaction systems are handled correctly"""

    reactions = chemkin('tests/test_xml/rxns_single.xml')
    expect = np.array([ -1.e4, -1.e4, 1.e4, 1.e4])
    assert all(reactions.reaction_rate([1,1,1,1],1000).astype(int) == expect.astype(int))
    return reactions


def test_reaction_rate():
    """
    test reaction rate claculation
    """
    reactions = test_from_xml()
    expect = np.array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
        -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    assert all(reactions.reaction_rate([1, 1, 1, 1, 1, 1],1000).astype(int) == expect.astype(int))

def test_construct_nonelementary():
    """
    test that chemkin object handles non-elementary reactions as expected
    """

    try:
        reactions = chemkin('tests/test_xml/rxns_pseudononelementary.xml')
        reactions.reaction_rate_T(np.ones(1), 1000)
    except NotImplementedError as e:
        assert type(e)==NotImplementedError
        print(e)

def test_invalid_reversible():
    """
    test that chemkin object handles invalid reversible input as expected
    """

    try:
        reactions = chemkin('tests/test_xml/rxns_invalid_reversible.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_repr():
    """
    test override of __repr__ function
    """
    reactions = test_from_xml()
    assert repr(reactions)=="chemkin('tests/test_xml/rxns.xml')"

def test_len():
    """
    test __len__ function
    """
    reactions = test_from_xml()
    assert len(reactions)==3

def test_str():
    """
    test __str__ function
    """
    reactions = test_from_xml()
    print(str(reactions))

def test_dimension_error():
    """
    check that inconsistencies in input params are handled correctly
    """
    try:
        reactions = chemkin('tests/test_xml/rxns_dimensionerror.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_bad_x():
    """
    test that errors with the shape of the concentration vector are handled correctly
    """
    reactions = test_from_xml()
    try:
        reactions.progress_rate([[1],[1]], 1000)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        reactions.reaction_rate([[1],[1]], 10000)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    
    try:
        reactions.reaction_rate([[1],[-1],[1],[1],[1],[1]], 1000)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
        
def test_long():
    reactions = chemkin('tests/test_xml/rxnset_long.xml')
    x = np.array([2., 1., 0.5, 1., 1., 0., 0., 0.25])
    T = 1500
    print(reactions.reaction_rate(x, T))
    expect = np.array([ -7.10942556e+12,  -1.58164670e+13,   2.20952666e+13,
        -1.65522031e+12,   1.12504921e+13,   0.00000000e+00,
         1.66470929e+13,  -2.54117388e+13])
    assert str(reactions.reaction_rate(x, T)) == str(expect)

def test_reversible():
    reactions = chemkin('tests/test_xml/rxns_reversible.xml')
    x = np.array([2., 1., 0.5, 1., 1., 0., 0., 0.25])
    T = 1500
    expect = np.array([  3.77689350e+14,  -3.81148508e+14,  -4.18872924e+14,   6.96275658e+12,
   3.07172197e+13,   3.93416746e+14,   1.66470991e+13,  -2.54117388e+13])
    assert str(reactions.reaction_rate(x, T)) == str(expect)

def test_singlereversible():
    reactions = chemkin('tests/test_xml/rxns_singlereversible.xml')
    x = np.array([2., 1., 0.5, 1.])
    T = 1500
    expect = np.array([  2.97178272e+14,   2.97178272e+14,  -2.97178272e+14,  -2.97178272e+14])
    assert str(reactions.reaction_rate(x, T)) == str(expect)

def test_mixedreversible():
    reactions = chemkin('tests/test_xml/rxnset_mixedreversible_long.xml')
    x = np.array([2., 1., 0.5, 1., 1., 0., 0., 0.25])
    T = 1500
    expect = np.array([ -1.51420660e+13,   1.22682375e+13,  -2.60415020e+13,   6.37742014e+12,
   3.13025562e+13,   0.00000000e+00,  1.66470929e+13,  -2.54117388e+13])
    assert str(reactions.reaction_rate(x, T)) == str(expect)
