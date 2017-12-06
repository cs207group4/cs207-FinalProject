import sys
import numpy as np

from pychemkin import InputParser


def test_correct():
    """test of standard usage"""
    input_ = InputParser('tests/test_xml/rxns.xml')
    assert len(input_) == 3
    assert len(eval(repr(input_))) == 3
    print(repr(input_))

def test_singlerxn():
    """test that input parser correctly handles the single reaction case"""
    input_ = InputParser('tests/test_xml/rxns_single.xml')
    assert len(input_) == 1
    assert np.array_equal(input_.nu_prod, np.array([[0],[0], [1], [1]]))
    assert np.array_equal(input_.nu_react, np.array([[1],[1], [0], [0]]))
    assert input_.rate_coeff_params[0]['k'] == 10000.0
    assert input_.rate_coeff_params[0]['type'] == 'Constant'
    print(repr(input_))

def test_extradata():
    """test that input parser correctly identifies too many reaction systems present"""
    try:
        input_ = InputParser('tests/test_xml/rxns_doubleddata.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_missingspecies():
    """test that input parser correctly identifies that species array is missing"""
    try:
        input_ = InputParser('tests/test_xml/rxns_missingspecies.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)

def test_missing_equation():
    """test that input parser correctly identifies that species array is missing"""
    try:
        input_ = InputParser('tests/test_xml/rxns_missing_attribute.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)

def test_bad_coeff():
    """test that input parser correctly identifies that species array is missing"""
    try:
        input_ = InputParser('tests/test_xml/rxns_bad_coeff.xml')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_missingcoeffdata():
    """test that input parser correctly identifies that coefficient data is missing"""
    try:
        input_ = InputParser('tests/test_xml/rxns_missingcoeffdata.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)
    try:
        input_ = InputParser('tests/test_xml/rxns_missing_k.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)
    try:
        input_ = InputParser('tests/test_xml/rxns_missing_A.xml')
    except AttributeError as e:
        assert type(e) == AttributeError
        print(e)
