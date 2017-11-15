import sys
import numpy as np
from pychemkin import InputParser, SQLParser, BackwardCoeffs, ReactionCoeffs

def test_CPcalc_single():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    Cp_r = bc.Cp_over_R(500.)
    
def test_Hcalc_single():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    Cp_r = bc.H_over_RT(500.)

def test_coeffs_single():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    rc_list = [ReactionCoeffs(**params) for params in ip.rate_coeff_params]
    for rc in rc_list:
        rc.set_params(T = 500)
    kf = np.array([rc.k_forward() for rc in rc_list])
    kb = bc.backward_coeffs(kf, 500)


def test_coeffs_None():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    bc.coeffs = None
    bc.H_over_RT(500)
    bc.coeffs = None
    bc.S_over_R(500)
    
def test_no_reversible():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    bc.backward_coeffs(np.array([]), 500)
    
def test_dimension():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    try:
        bc.backward_coeffs(np.array([1, 2, 3]), 500)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)