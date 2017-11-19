import sys
import numpy as np
from pychemkin import InputParser, SQLParser, BackwardCoeffs, ReactionCoeffs

def test_CPcalc_single():
    """
    Test that single reaction case is handled correctly
    """
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser()
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species)
    Cp_r = bc.Cp_over_R(500.)
    assert np.all(np.isclose(Cp_r,np.array([ 2.5,3.73848592,2.55540662,3.54595678])))

def test_Hcalc_single():
    """
    Test that single reaction case is handled correctly
    """
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser()
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species)
    H_ = bc.H_over_RT(500.)
    assert np.all(np.isclose(H_,np.array([ 53.4473198,1.46400033,60.98145063, 10.90244139])))

def test_coeffs_single():
    """
    Test that single reaction case is handled correctly
    """
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser()
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species)
    rc_list = [ReactionCoeffs(**params) for params in ip.rate_coeff_params]
    for rc in rc_list:
        rc.set_params(T = 500)
    kf = np.array([rc.k_forward() for rc in rc_list])
    kb = bc.backward_coeffs(kf, 500)
    assert np.isclose(kb[0],6.03552730e+18)


def test_coeffs_None():
    """
    Check that coefficients are being set correctly
    """
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser()
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species)
    bc.coeffs = None
    H_ = bc.H_over_RT(500)
    assert np.all(np.isclose(H_,np.array([ 53.4473198,1.46400033,60.98145063, 10.90244139])))
    bc.coeffs = None
    S_ = bc.S_over_R(500)
    assert np.all(np.isclose(S_, np.array([ 15.08983739,26.54400129,20.70986963,23.94106991])))

def test_no_reversible():
    """
    Check that the empty case can be handled
    """
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser()
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species)
    assert len(bc.backward_coeffs(np.array([]), 500))==0

def test_dimension():
    """
    Check that incompatible input array sizes are correctly caught
    """
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser()
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species)
    try:
        bc.backward_coeffs(np.array([1, 2, 3]), 500)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
