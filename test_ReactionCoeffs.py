from chemkin import *
def test_str():
    rc = ReactionCoeffs('Constant', k = 1e3)
    correct = "Constant Reaction Coeffs:"
    assert str(rc)[:len(correct)] == correct
def test_repr():
    rc = ReactionCoeffs('modifiedArrhenius', A = 1e7, E=1e3, T=1e2, b=0.5)
    rc1 = eval(repr(rc))
    assert rc == rc1
def test_insufficient():
    rc = ReactionCoeffs('modifiedArrhenius', A = 1e7, E=1e3, b=0.5)
    try:
        rc.kval()
    except ValueError as e:
        assert type(e) == ValueError