from chemkin import *
def test_from_xml():
    """
    testing normal reaction rate calculations from xml input
    """
    reactions = chemkin.from_xml('rxns.xml')
    reactions.set_rc_params(T=1000)
    expect = np.array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
        -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    assert all(reactions.reaction_rate([1, 1, 1, 1, 1, 1]).astype(int) == expect.astype(int))
    return reactions


def test_from_singlerxnxml():

    """testing that single reaction systems are handled correctly"""

    reactions = chemkin.from_xml('test_xml/rxns_single.xml')
    reactions.set_rc_params(T=1000)
    expect = np.array([ -1.e4, -1.e4, 1.e4, 1.e4])
    assert all(reactions.reaction_rate([1,1,1,1]).astype(int) == expect.astype(int))
    return reactions


def test_reaction_rate_with_T():
    """
    test wrapper function for specifying temp and concentration in the same function
    """
    reactions = test_from_xml()
    expect = np.array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
        -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    assert all(reactions.reaction_rate_T([1, 1, 1, 1, 1, 1],1000).astype(int) == expect.astype(int))

def test_construct():
    """
    test direct input of reaction system params into instantiation of chemkin object
    """
    rc_list = [
        ReactionCoeffs(type = "modifiedArrhenius", A=1e8, b=0.5, E=5e4),
        ReactionCoeffs(type = "Constant", k=1e4),
        ReactionCoeffs(type = "Arrhenius", A=1e7, E=1e4),
    ]
    v1 = np.array([
        [2,0,0],
        [1,0,1],
        [0,1,0],
        [0,1,0],
        [0,0,1]
    ])
    v2 = np.array([
        [1,0,0],
        [0,1,0],
        [2,0,1],
        [0,0,1],
        [0,1,0]
    ])
    x = np.array([
        2,1,0.5,1,1
    ]).reshape(-1,1)
    expect = dict()
    expect[750] = np.array([-3607077.87280406, -5613545.18362079,  9220623.05642485,
        2006467.31081673, -2006467.31081673])
    expect[1500] = np.array([ -2.81117621e+08,  -2.85597559e+08,   5.66715180e+08,
         4.47993847e+06,  -4.47993847e+06])
    expect[2500] = np.array([ -1.80426143e+09,  -1.81043736e+09,   3.61469878e+09,
         6.17593098e+06,  -6.17593098e+06])
    rxn_types = ['Elementary']*5
    reversible = ['no']*5
    reactions = chemkin(v1,v2,rc_list,rxn_types, reversible)

    for t in [750,1500,2500]:
        reactions.set_rc_params(T=t)
        assert str(reactions.reaction_rate(x)) == str(expect[t])
    reactions = chemkin.init_const_rc(v1,v2,[1,1,1],['Elementary']*5, ['no']*5, ['A','B','X','Y','Z'])
    print(reactions)

def test_construct_nonelementary():
    """
    test that chemkin object handles non-elementary reactions as expected
    """
    #these are just dummy values. They would not be valid entries for non-elementary reactions
    rc_list = [
        ReactionCoeffs(type = "modifiedArrhenius", A=1e8, b=0.5, E=5e4),
        ReactionCoeffs(type = "Constant", k=1e4),
        ReactionCoeffs(type = "Arrhenius", A=1e7, E=1e4),
    ]
    v1 = np.array([
        [2,0,0],
        [1,0,1],
        [0,1,0],
        [0,1,0],
        [0,0,1]
    ])
    v2 = np.array([
        [1,0,0],
        [0,1,0],
        [2,0,1],
        [0,0,1],
        [0,1,0]
    ])
    x = np.array([
        2,1,0.5,1,1
    ]).reshape(-1,1)
    expect = dict()
    expect[750] = np.array([-3607077.87280406, -5613545.18362079,  9220623.05642485,
        2006467.31081673, -2006467.31081673])
    expect[1500] = np.array([ -2.81117621e+08,  -2.85597559e+08,   5.66715180e+08,
         4.47993847e+06,  -4.47993847e+06])
    expect[2500] = np.array([ -1.80426143e+09,  -1.81043736e+09,   3.61469878e+09,
         6.17593098e+06,  -6.17593098e+06])
    rxn_types = ['Elementary', 'Hello Kitty', 'Elementary', 'Elementary', 'Elementary']
    reversible = ['no']*5
    reactions = chemkin(v1,v2,rc_list,rxn_types, reversible)
    try:
        reactions.reaction_rate_T(np.ones(5), 1000)
    except NotImplementedError as e:
        assert type(e)==NotImplementedError
        print(e)

def test_str_repr():
    """
    test override of __str__ and __repr__ functions
    """
    reactions = test_from_xml()
    reactions_1 = eval(repr(reactions))
    print(reactions_1)
    print(reactions)
    assert 3==len(reactions_1)
    reactions.species=None
    print(reactions)

def test_dimension_error():
    """
    check that inconsistencies in input params are handled correctly
    """
    try:
        reactions = chemkin.init_const_rc([[1,2],[1,2]],[[1,2],[1,2]],[1,1,1],['Elementary']*5,['no']*5,['A','B','X','Y','Z'])
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_bad_x():
    """
    test that errors with the shape of the concentration vector are handled correctly
    """
    reactions = test_from_xml()
    try:
        reactions.progress_rate([[1],[1]])
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        reactions.reaction_rate([[1],[1]])
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        reactions.reaction_rate([[1],[-1],[1],[1],[1],[1]])
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        reactions.reaction_rate_T([[1],[-1],[1],[1],[1],[1]],-1)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
