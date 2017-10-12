from chemkin import *
def test_str():
    rc = reaction_coeffs('Constant', k = 1e3)
    correct = "Constant Reaction Coeffs:"
    assert str(rc)[:len(correct)] == correct
def test_repr():
    rc = reaction_coeffs('modifiedArrhenius', A = 1e7, E=1e3, T=1e2, b=0.5)
    rc1 = eval(repr(rc))
    assert rc == rc1