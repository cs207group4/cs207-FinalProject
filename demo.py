import numpy as np
import chemkin

rxns = chemkin.chemkin.from_xml('rxnset_long.xml')
x = np.array([2., 1., 0.5, 1., 1., 0., 0., 0.25])
T = 1500
print(rxns.reaction_rate_T(x, T))
