from pychemkin import ReactionCoeffs


def test_str():
    """
    testing override of __str__ for ReactionCoeffs class
    """
    rc = ReactionCoeffs('Constant', k = 1e3)
    correct = "Constant Reaction Coeffs:"
    assert str(rc)[:len(correct)] == correct

def test_repr():
    """
    testing override of __repr__ for ReactionCoeffs class
    """
    rc = ReactionCoeffs('modifiedArrhenius', A = 1e7, E=1e3, T=1e2, b=0.5)
    rc1 = eval(repr(rc))
    assert rc == rc1

def test_insufficient():
    """
    test when reaction coefficient params are not suitable for reaction type
    """
    rc = ReactionCoeffs('modifiedArrhenius', A = 1e7, E=1e3, b=0.5)
    try:
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
    rc = ReactionCoeffs('Arrhenius', b=0.5)
    try:
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
    rc = ReactionCoeffs('Constant', A = 1e7, E=1e3, b=0.5)
    try:
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_bad_arg():
    """
    test for incorrect or extraneous arguments
    """
    rc = ReactionCoeffs('modifiedArrhenius', A = 1e7, E=1e3, b=(0.5,8j), T=1000)
    try:
        rc.k_forward()
    except TypeError as e:
        assert type(e) == TypeError
        print(e)
    try:
        rc.set_params(b=0.5,T=-1)
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        rc.set_params(b=0.5,T=1,A=-1)
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    rc = ReactionCoeffs('Arrhenius', A = 1e7, E=1e3, T=-1000)
    try:
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        rc.set_params(A=-1,T=1)
        rc.k_forward()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_not_implement():
    """test of hooks (things not yet implemented)"""
    try:
        rc = ReactionCoeffs('HelloWorld', A = 1e7, E=1e3, b=0.5)
    except NotImplementedError as e:
        assert type(e) == NotImplementedError
        print(e)
