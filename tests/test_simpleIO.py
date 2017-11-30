from pychemkin import chemkin, ChemSolver, simpleIO
import numpy as np

def test_simpleIO():
    chem2 = chemkin("tests/test_xml/rxns_reversible.xml")
    # initial concentration
    x_init = np.ones(8)
    T = 2000
    t_final = 3.e-13
    t_step = 1.e-15
    cs = ChemSolver(chem2).solve(x_init, T, t_final, t_step, algorithm = 'lsoda')
    r1 = cs.get_results()

    simpleIO('tests/test_data/test.pkl').to_pickle(cs)
    cs2 = simpleIO('tests/test_data/test.pkl').read_pickle()
    r2 = cs2.get_results()

    for i in range(len(r1)):
        assert np.any(r1[i] == r2[i])

    assert repr(cs.chem) == repr(cs2.chem)
    
def test_wrong_file_extension():
    try:
        simpleIO('tests/test_data/test.csv')
    except ValueError as e:
        assert type(e) == ValueError