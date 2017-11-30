'''

pychemkin
~~~~~~~~~

A python library for computing chemical kinetics. For details, please refer to the user manual. The main class user should use is chemkin.

:copyright: (c) 2017 by cs207 group4
:license: GPLv3, see LICENSE.txt for details

'''

from .simpleIO import simpleIO
from .InputParser import InputParser
from .SQLParser import SQLParser
from .BackwardCoeffs import BackwardCoeffs
from .ReactionCoeffs import ReactionCoeffs
from .chemkin import chemkin
from .ChemSolver import ChemSolver
from .ChemViz import ChemViz