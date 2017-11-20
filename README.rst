Python Chemical Kinetics Library
==================================
pychemkin is a Python 3 library that computes the reaction rates of all species participating
in a system of elementary reactions.

The package can solve for the reaction rates of a system of elementary reactions. The
number of reactions and species is arbitrary. For each system of reactions, the user supplies
the species participating in the reactions, the chemical equations, the stoichiometric
coefficients for the reactants and products, and the rate coefficient parameters (e.g., E and
A for Arrhenius rates). For a given system, the user can then specify a temperature and
a vector of species concentrations in order to return the reaction rates in the form of a
NumPy array. Rate coefficients and reaction progress rates can also be retrieved.

Further documentation is in 'docs/pychemkin.pdf' within the source file.
