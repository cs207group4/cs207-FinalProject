
import sys
import numpy as np
from pychemkin import *

def test_CPcalc_single():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    Cp_r = bc.Cp_over_R(500.)
    
def test_Hcalc_single():
    pass

def test_coeffs_single():
    ip = InputParser('tests/test_xml/rxns_singlereversible.xml')
    sql = SQLParser('pychemkin/data/thermo.sqlite')
    bc = BackwardCoeffs(ip.nu_react, ip.nu_prod, ip.species, sql)
    rc_list = [ReactionCoeffs(**params) for params in ip.rate_coeff_params]
    for rc in rc_list:
        rc.set_params(T = 500)
    kf = np.array([rc.k_forward() for rc in rc_list])
    kb = bc.backward_coeffs(kf, 500)


def test_coeffs():
    pass
