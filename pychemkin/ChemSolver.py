import numpy as np
import pandas as pd
import warnings
from scipy.integrate import ode
from .chemkin import chemkin


class ChemSolver:
    '''
    The ChemSolver module wraps around SciPy's ODE solver to obtain the evolution of concentrations and reaction
    rates over time

    METHODS and ATTRIBUTES
    ========
    After initialization, user could call:
     - solve(y0, T, t_span, t_eval=None, **options): method to solve the differential equations given an initial
       value of species concentrations:
           dy / dt = chemkin.reaction_rate(y, T)
           y(t0) = y0
     - get_results(return_reaction_rate=True): method to return the solution of ODEs, which includes an array of
       time points, an array of solution values (species concentrations) at each time point, and an array of
       reaction rates (optional)
     - save_results(file_name): method to save the solution of ODEs to a csv or hdf5 file.
     - load_results(file_name): method to load the solution of ODEs from a csv or hdf5 file storing the data
     - grid_search(y0s, Ts, t_span, return_reaction_rate=True, **options): method to solve ODE at different
       combinations of starting concentrations and temperatures
     - get_grid_search_result(): method to return the grid search starting conditions and results

    EXAMPLES
    =========
    >>> chem = chemkin("tests/test_xml/rxns.xml")
    Finished reading xml input file
    >>> y0 = np.ones(len(chem.species))
    >>> T = 300
    >>> t1 = 0.001
    >>> dt = 2e-4
    >>> cs = ChemSolver(chem).solve(y0, T, t1, dt, method='lsoda')
    >>> t, y, rr = cs.get_results()
    >>> y[0]
    array([ 1.        ,  1.16784753,  1.28789036,  1.3780344 ,  1.44823276])
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
        self.finished = False

    def __dy_dt(self, t, y):
        r = self.chem._progress_rate_default_T(y.reshape(-1,1),self.kf,self.kb)
        return np.sum(r * self.nu_diff_matrix, axis=1)

    def solve(self, y0, T, t1, dt, algorithm='lsoda', **options):
        '''Solve ODEs

        INPUTS
        ======
        y0: array_like, shape(n,), required
            A length-n vector specifying initial concentrations of the n species
        T: float, required
            Temperature in K
        t1: float, required
            Integration end point. The solver starts with t=0 and integrates until it reaches t=t1.
        dt: float, required
            Integration step.
        algorithm: string, optional
            Integration algorithm. default value is 'lsoda'
        **options:
            Options passed to scipy.integrate.ode.set_integrator
        '''
        self.T = T
        # get variables not changed in ODE
        _,self.kf, self.kb = self.chem._progress_rate_init(y0, T)
        self.nu_diff_matrix = self.chem.nu_prod-self.chem.nu_react
        r = ode(self.__dy_dt).set_integrator(algorithm, **options)
        r.set_initial_value(y0, 0)
        self._t = [0]
        self._y = [y0]
        total_pts = t1/dt
        while r.successful() and len(self._t)<total_pts:
            self._t.append(r.t + dt)
            self._y.append(r.integrate(r.t + dt))

        #check whether integration completed properly
        if len(self._t) == total_pts:
            self.finished = True
        else:
            print("WARNING: Solver appears to have stopped prematurely. Try increasing nsteps or decreasing dt.")
        self._t = np.array(self._t)
        self._y = np.array(self._y).transpose()
        self._sol = True
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
        t = self._t
        y = self._y
        if return_reaction_rate and self.reaction_rate is None:
            self.reaction_rate = np.concatenate([self.__dy_dt(0, y[:, i]).\
                                                 reshape((y.shape[0], 1)) for i in range(y.shape[1])], axis=1)
        if return_reaction_rate:
            return t, y, self.reaction_rate
        else:
            return t, y, None

    def is_equilibrium(self,tol=1.):
        '''If the reaction system has reached equilibrium

        INPUT
        =======
        tol: allowed error. If all derivatives less than tol, system will be considered equilibrium

        RETURN
        =======
        is_equilibrium: bool. Indicate equilibrium or not
        '''

        if self._sol is None:
            raise ValueError('ODE not solved.')
        return np.all(np.abs(self.__dy_dt(0, self._y[:,-1]))<tol)

    def to_df(self):
        '''Save the solution of ODEs to a pandas dataframe

        RETURNS
        =======
        df: pandas data frame containing the solution

        '''
        if self._sol is None:
            raise ValueError('ODE not solved.')
        t = self._t
        y = self._y
        if self.reaction_rate is None:
            self.reaction_rate = np.concatenate([self.__dy_dt(0, y[:, i]).\
                                                 reshape((y.shape[0], 1)) for i in range(y.shape[1])], axis=1)
        rr = self.reaction_rate
        cols = ['t'] + ['{}-Concentration'.format(s) for s in self.chem.species] \
        + ['{}-Reaction_rate'.format(s) for s in self.chem.species] + ['T']
        values = np.concatenate((t.reshape(1, len(t)), y.reshape(len(y), len(t)), \
                                 rr.reshape(len(rr), len(t)), self.T * np.ones((1, (len(t))))), axis=0)
        return pd.DataFrame.from_records(values, cols).transpose()

    def save_results(self, file_name):
        '''Save the solution of ODEs to a csv or hdf5 file

        INPUT
        ======
        file_name: string, required
            csv or hdf5 file to save the solution of ODEs
            The file extension should be either .csv or .h5
        '''
        file_name = file_name.strip()
        if not ('.csv' == file_name[-4:] or '.h5' == file_name[-3:]):
            raise ValueError('Only csv and hdf5 are supported.')
        df = self.to_df()
        if '.csv' == file_name[-4:]:
            df.to_csv(file_name, index=False)
        else:
            df.to_hdf(file_name, 'df',index=False)

    def load_results(self, file_name):
        '''Load the solution of ODEs from a csv or hdf5 file

        INPUT
        ======
        file_name: string, required
            csv or hdf5 file containing the ODE solutions
            The file extension should be either .csv or .h5
        '''
        file_name = file_name.strip()
        if not ('.csv' == file_name[-4:] or '.h5' == file_name[-3:]):
            raise ValueError('Only csv and hdf5 are supported.')
        if '.csv' == file_name[-4:]:
            df = pd.read_csv(file_name)
        else:
            df = pd.read_hdf(file_name, 'df')
        self._sol = True
        self._t = df.t.values
        self._y = df.iloc[:, 1:1+len(self.chem.species)].values.transpose()
        self.reaction_rate = df.iloc[:, 1+len(self.chem.species):1+2*len(self.chem.species)].values.transpose()
        self.T = df['T'][0]
        _,self.kf, self.kb = self.chem._progress_rate_init(self._y.transpose()[0], self.T)
        self.nu_diff_matrix = self.chem.nu_prod-self.chem.nu_react
        return self

    def grid_search(self, y0s, Ts, t1, dt, method='lsoda', return_reaction_rate=True, **options):
        '''Solve ODEs at user specified combinations of starting concentrations and temperatures

        INPUTS
        =======
        y0s: array_like, shape(n_y_cond, n), required
            A n_y_cond X n array specifying n_y_cond groups of initial concentrations of the n species
        T: array_like, shape(n_T_cond,), required
            A length-n_T_cond vector specifying n_T_cond temperatures in K
        t1: float, required
            Integration end point. The solver starts with t=0 and integrates until it reaches t=t1.
        dt: float, required
            Integration step.
        method: string, optional
            Integration algorithm. default value is 'lsoda'
        return_reaction_rate: bool, optional
            If True (default), calculate reaction rates
        **options:
            Options passed to scipy.integrate.set_integrator
        '''
        self.grid_result = {T:[self.solve(y0, T, t1, dt, method, **options).\
                               get_results(return_reaction_rate) for y0 in y0s] for T in Ts}
        self.grid_condition = [Ts, y0s]
        self._sol = None
        self.reaction_rate = None
        return self

    def get_grid_results(self):
        '''Return the grid search starting conditions and results
        '''
        return self.grid_condition, self.grid_result
