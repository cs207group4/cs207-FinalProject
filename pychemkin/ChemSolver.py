import numpy as np
from scipy.integrate import solve_ivp
from .chemkin import chemkin

class ChemSolver:
    def __init__(self, chem):
        self.chem = chem
        self._sol = None
        self.reaction_rate = None
        self.T = None
        self.grid_condition = None
        self.grid_result = None
    
    def solve(self, y0, T, t_span, t_eval=None, **options):
        self.T = T
        def _ode(t, y0):
            return self.chem.reaction_rate(y0, T)
        self._sol = solve_ivp(_ode, t_span=t_span, y0=y0, t_eval=t_eval, **options)
        self.reaction_rate = None
        return self
    
    def get_results(self, return_reaction_rate=True):
        if self._sol is None:
            raise ValueError('ODE not solved.')
        t = self._sol.t
        y = self._sol.y
        if return_reaction_rate and self.reaction_rate is None:
            self.reaction_rate = np.concatenate([self.chem.reaction_rate(y[:, i], self.T).\
                                                 reshape((y.shape[0], 1)) for i in range(y.shape[1])], axis=1)
        if return_reaction_rate:
            return t, y, self.reaction_rate
        else:
            return t, y, None
        
    def grid_search(self, y0s, Ts, t_span, return_reaction_rate=True, **options):
        self.grid_result = {T:[self.solve(y0, T, t_span, **options).get_results(return_reaction_rate) \
                            for y0 in y0s] for T in Ts}
        self.grid_condition = [Ts, y0s] 
        self._sol = None
        self.reaction_rate = None
        return self
    
    def get_grid_results(self):
        return self.grid_condition, self.grid_result