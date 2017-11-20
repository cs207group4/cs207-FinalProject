from pychemkin import chemkin, ChemSolver
import numpy as np

def test_solve():
    chem = chemkin("tests/test_xml/rxns.xml")
    y0 = np.ones(len(chem.species))
    T = 300
    t_span = (0, 0.003)
    cs = ChemSolver(chem).solve(y0, T, t_span, method='RK23')
    t, y, rr = cs.get_results(return_reaction_rate=False)
    assert rr is None
    
    t, y, rr = cs.get_results()
    
    expect_y_s0 = np.array([ 1., 1.15629864, 1.29624502, 1.41687655, 1.5192147, 1.60512333, \
                            1.67676943, 1.73628759, 1.75927183])
    
    assert rr.shape == (6, 9)
    assert str(y[0]) == str(expect_y_s0)
    print(y)
    
def test_grid_search():
    chem = chemkin("tests/test_xml/rxns.xml")
    cs = ChemSolver(chem)
    t_span = (0, 0.003)
    y0s = [np.ones(len(chem.species))*i for i in range(1, 4)]
    Ts = [300, 400]
    gs = cs.grid_search(y0s, Ts, t_span, method='RK23')
    _y1 = gs.get_grid_results()[1][300][0][1]
    _y2 = cs.solve(y0s[0], Ts[0], t_span, method='RK23').get_results()[1]
    assert str(_y1) == str(_y2)
    print(_y1)
    
def test_OED_not_solved():
    chem = chemkin("tests/test_xml/rxns.xml")
    cs = ChemSolver(chem)
    try:
        t, y, rr = cs.get_results()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)