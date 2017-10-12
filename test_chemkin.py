from chemkin import *
def test_from_xml():
    reactions = chemkin.from_xml('rxns.xml')
    reactions.set_rc_params(T=1000)
    expect = np.array([ -6.28889929e+06,   6.28989929e+06,   6.82761528e+06,
        -2.70357993e+05,   1.00000000e+03,  -6.55925729e+06])
    assert all(reactions.reaction_rate([[1],[1],[1],[1],[1],[1]]).astype(int) == expect.astype(int))

def test_construct():
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
    reactions = chemkin(v1,v2,rc_list)
    for t in [750,1500,2500]:
        reactions.set_rc_params(T=t)
        assert str(reactions.reaction_rate(x)) == str(expect[t])