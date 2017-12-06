from pychemkin import chemkin, ChemSolver,ChemViz
import numpy as np
import os
def test_html_irrev():
    chem = chemkin("tests/test_xml/rxns.xml")
    y0 = np.ones(len(chem.species))
    T = 300
    t1 = 0.003
    dt = 5e-4
    cs = ChemSolver(chem).solve(y0, T, t1, dt, method='lsoda')
    t, y, rr = cs.get_results(return_reaction_rate=False)
    ChemViz(cs).html_report('report1.html')
    assert os.path.isfile('report1.html')

def test_plot_gs():
    chem = chemkin("tests/test_xml/rxns.xml")
    t1 = 0.003
    dt = 5e-4
    y0s = [np.ones(len(chem.species))*i for i in range(1, 4)]
    Ts = [300, 400]
    gs = ChemSolver(chem).grid_solve(y0s, Ts, t1, dt, algorithm='lsoda')
    ChemViz(gs).plot_gridtime_series('reactionrate')
    ChemViz(gs).plot_gridtime_series('concentration')
    try:
        ChemViz(gs).plot_gridtime_series('reactionrate',tmin=-1)
    except ValueError as e:
        assert type(e) == ValueError
        print(e)


def test_reversible_report():
    x_init = np.ones(8)
    T = 1000

    # integration end time
    t_max = 5.e-13

    # step size
    dt = 1.e-16

    cs = ChemSolver(chemkin('tests/test_xml/rxns_reversible.xml'))
    try:
        ChemViz(cs).html_report('report2.html')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)
    cs.solve(x_init, T, t_max, dt)
    ChemViz(cs).html_report('report2.html')
    assert os.path.isfile('report2.html')

def test_demo_report():
    x_init = np.ones(8)
    T = 1000

    # integration end time
    t_max = 5.e-13

    # step size
    dt = 1.e-16

    cs = ChemSolver(chemkin('tests/test_xml/demo.xml')).solve(x_init, T, t_max, dt)
    ChemViz(cs).html_report('report3.html')
    assert os.path.isfile('report3.html')

def test_no_grid():
    x_init = np.ones(8)
    T = 1000

    # integration end time
    t_max = 5.e-13

    # step size
    dt = 1.e-16

    cs = ChemSolver(chemkin('tests/test_xml/demo.xml')).solve(x_init, T, t_max, dt)
    try:
        ChemViz(cs).plot_gridtime_series('reactionrate')
    except ValueError as e:
        assert type(e) == ValueError
        print(e)