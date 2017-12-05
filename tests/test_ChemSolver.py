from pychemkin import chemkin, ChemSolver
import numpy as np

def test_solve():
    chem = chemkin("tests/test_xml/rxns.xml")
    y0 = np.ones(len(chem.species))
    T = 300
    t1 = 0.003
    dt = 5e-4
    cs = ChemSolver(chem).solve(y0, T, t1, dt, algorithm='lsoda')
    t, y, rr = cs.get_results(return_reaction_rate=False)
    assert rr is None
    t, y, rr = cs.get_results()
    print(y)

def test_grid_search():
    chem = chemkin("tests/test_xml/rxns.xml")
    t1 = 0.003
    dt = 5e-4
    y0s = [np.ones(len(chem.species))*i for i in range(1, 4)]
    Ts = [300, 400]
    gs = ChemSolver(chem).grid_solve(y0s, Ts, t1, dt, algorithm='lsoda')
    _y1 = gs.get_grid_results()[1][300][0][1]
    _y2 = ChemSolver(chem).solve(y0s[0], Ts[0], t1, dt, algorithm='lsoda').get_results()[1]
    assert str(_y1) == str(_y2)
    print(_y1)

def test_ODE_not_solved():
    chem = chemkin("tests/test_xml/rxns.xml")
    cs = ChemSolver(chem)
    try:
        t, y, rr = cs.get_results()
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        cs.save_results('tests/test_data/test.csv')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)

def test_IO():
    chem = chemkin("tests/test_xml/rxns.xml")
    y0 = np.ones(len(chem.species))
    T = 300
    t1 = 0.003
    dt = 5e-4
    cs = ChemSolver(chem).solve(y0, T, t1, dt, algorithm='lsoda')
    cs.save_results('tests/test_data/test.csv')
    cs.save_results('tests/test_data/test.h5')
    cs1 = ChemSolver(chem).load_results('tests/test_data/test.csv')
    cs2 = ChemSolver(chem).load_results('tests/test_data/test.h5')
    assert str(cs.get_results()) == str(cs1.get_results())
    assert str(cs1.get_results()) == str(cs2.get_results())

def test_wrong_file_name():
    chem = chemkin("tests/test_xml/rxns.xml")
    cs = ChemSolver(chem)
    try:
        cs.save_results('test.pdf')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    try:
        cs.load_results('test.pdf')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
