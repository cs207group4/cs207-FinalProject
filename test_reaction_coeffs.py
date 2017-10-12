from reaction_coeffs import *
def test_str():
    rc = reaction_coeffs('Constant', k = 1e3)
    assert str(rc) == "Constant Reaction Coeffs: {'k': 1000.0, 'R': 8.314}"
def test_repr():
    rc = reaction_coeffs('modifiedArrhenius', A = 1e7, E=1e3, T=1e2, b=0.5)
    rc1 = eval(repr(rc))
    assert rc == rc1