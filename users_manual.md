# User Manual
## Introduction: 

`chemkin` is a Python library that computes the reaction rates of all species participating in a system of elementary, irreversible reactions. 

### Key chemical concepts and terminology

A system consisting of $M$ elementary reactions involving $N$ species has the general form 
\begin{align}
  \sum_{i=1}^{N}{\nu_{ij}^{\prime}\mathcal{S}_{i}} \longrightarrow 
  \sum_{i=1}^{N}{\nu_{ij}^{\prime\prime}\mathcal{S}_{i}}, \qquad j = 1, \ldots, M
\end{align}

$S_i$ is the $i$th specie in the system, $\nu_{ij}^{\prime}$ is its stoichiometric coefficient (dimensionless) on the reactants side of the $j$th reaction, and ${\nu_{ij}^{\prime\prime}$ is its stoichiometric coefficent (dimensionless) on the product side for the $j$th reaction. 


### Features

The package can solve for the reaction rates of a system with an arbitrary number of species and elementary reactions. The reaction rate coefficient $k$ for each reaction is assumed to take one of three possible forms:
1. $k=$ constant
2. Arrhenius: $k=A\exp(-\frac{E}{RT})$, where $A$ is the pre-factor, $E$ is the activation energy, $R$ is the universal gas constant, and $T$ is the temperature. 
3. Modified Arrhenius: $k=AT^b\exp(-\frac{E}{RT})$, where $A$ is the pre-factor, $E$ is the activation energy, $R$ is the universal gas constant, $T$ is the temperature, and $b$ is the temperature scaling parameter. 

#### Input
#### Output
## Installation: 

Describe where the code can be found and downloaded. Tell the user how to run the test suite. We are not releasing this code as a package yet, but when we do that this section will include instructions how how to install the package.

## Basic Usage and Examples: 

### Input format

Provide a few examples on using your software in some common situations. You may want to show how the code works with a small set of reactions.

