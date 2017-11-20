import numpy as np
import warnings
from scipy.integrate import solve_ivp
from .chemkin import chemkin

class ChemSolver:
    '''
    The ChemSolver module wrap around SciPy's ODE solver to obtain the evolution of concentrations and reaction 
    rates over time
    
    METHODS and ATTRIBUTES
    ========
    After initialization, user could call:
     - solve(y0, T, t_span, t_eval=None, **options): method to solve the differential equations given an initial value
       of species concentrations:
           dy / dt = chemkin.reaction_rate(y, T)
           y(t0) = y0
     - get_results(return_reaction_rate=True): method to return the solution of ODEs, which includes an array of time 
       points, an array of solution values (species concentrations) at each time point, and an array of species reaction 
       rates (optional)
     - grid_search(y0s, Ts, t_span, return_reaction_rate=True, **options): method to solve ODE at different combinations
       of starting concentrations and temperatures
     - get_grid_search_result(): method to return the grid search starting conditions and results
     
    EXAMPLES
    =========
    >>> chem = chemkin("tests/test_xml/rxns.xml")
    Finished reading xml input file
    >>> y0 = np.ones(len(chem.species))
    >>> T = 300
    >>> t_span = (0, 0.003)
    >>> cs = ChemSolver(chem).solve(y0, T, t_span, method='RK23')
    >>> t, y, rr = cs.get_results()
    >>> y[0]
    array([ 1.        ,  1.15629864,  1.29624502,  1.41687655,  1.5192147 ,
            1.60512333,  1.67676943,  1.73628759,  1.75927183])
    '''
    def __init__(self, chem):
        '''
        INPUT
        =====
        chem: chemkin object, required
        '''
        self.chem = chem
        self._sol = None
        self.reaction_rate = None
        self.T = None
        self.grid_condition = None
        self.grid_result = None
    
    def solve(self, y0, T, t_span, t_eval=None, **options):
        '''Solve ODEs
        
        INPUTS
        ======
        y0: array_like, shape(n,), required
            A length-n vector specifying initial concentrations of the n species
        T: float, required
            Temperature in K
        t_span: 2-tuple of floats, required
            Interval of integration (t0, tf). The solver starts with t=t0 and integrates until it reaches t=tf.
        t_eval: array_like or None, optional
            Times at which to store the computed solution, must be sorted and lie within t_span. 
            If None (default), use points selected by a solver.
        **options:
            Options passed to scipy.integrate.solve_ivp
            see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html for details
        '''
        self.T = T
        def _ode(t, y0):
            return self.chem.reaction_rate(y0, T)
        self._sol = solve_ivp(_ode, t_span=t_span, y0=y0, t_eval=t_eval, **options)
        if -1 == self._sol.status:
            warnings.warn('Integration step failed.', Warning)
        elif 1 == self._sol.status:
            warnings.warn('A termination event occurred.', Warning)
        self.reaction_rate = None
        return self
    
    def get_results(self, return_reaction_rate=True):
        '''Return ODE solutions
        
        INPUT
        =====
        return_reaction_rate: bool, optional
            If True (default), return reaction rates
        
        RETURNS
        =======
        t: ndarray, shape(n_points,)
            Time points.
        y: ndarray, shape(n, n_points)
            Solution values (concentrations) at t
        reaction_rate: ndarray, shape(n, n_points)
            Reaction rates at t
        '''
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
        '''Solve ODEs at user specified combinations of starting concentrations and temperatures
        
        INPUTS
        =======
        y0s: array_like, shape(n_y_cond, n), required
            A n_y_cond X n array specifying n_y_cond groups of initial concentrations of the n species
        T: array_like, shape(n_T_cond,), required
            A length-n_T_cond vector specifying n_T_cond temperatures in K
        t_span: 2-tuple of floats, required
            Interval of integration (t0, tf). The solver starts with t=t0 and integrates until it reaches t=tf.
        return_reaction_rate: bool, optional
            If True (default), calculate reaction rates
        **options:
            Options passed to scipy.integrate.solve_ivp
            see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html for details
        '''
        self.grid_result = {T:[self.solve(y0, T, t_span, **options).get_results(return_reaction_rate) \
                            for y0 in y0s] for T in Ts}
        self.grid_condition = [Ts, y0s] 
        self._sol = None
        self.reaction_rate = None
        return self
    
    def get_grid_results(self):
        '''Return the grid search starting conditions and results
        '''
        return self.grid_condition, self.grid_result