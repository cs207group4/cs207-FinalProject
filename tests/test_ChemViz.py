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
