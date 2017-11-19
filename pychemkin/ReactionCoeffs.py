import numpy as np

class ReactionCoeffs:
    """
    A class for computing forward reaction coefficients

    set_params: sets reaction rate coefficient parameters
    k_forward: compute forward reaction rate coefficient

    EXAMPLES
    =========
    >>> rc = ReactionCoeffs('Constant', k = 1e3)
    >>> rc.k_forward()
    1000.0
    >>> rc = ReactionCoeffs('Arrhenius', A = 1e7, E=1e3)
    >>> rc.set_params(T=1e2)
    >>> rc.k_forward()
    3003549.0889639612
    """

    def __init__(self, type, **kwargs):
        '''
        INPUT
        =====
        type: string, required
              Type of the reaction rate coefficient of interest (e.g. "Constant", "Arrhenius", etc)
        '''
        self.__rtype = str(type)
        self.__params = kwargs
        if 'R' not in self.__params:
            self.__params['R'] = 8.314

    def set_params(self, **kwargs):
        for param in kwargs:
            self.__params[param] = kwargs[param]

    def __str__(self):
        return self.__rtype + " Reaction Coeffs: " + str(self.__params)

    def __repr__(self):
        class_name = type(self).__name__
        params_str = ""
        for param in self.__params:
            params_str += " ,"+param+" = "+str(self.__params[param])
        return str(class_name)+'("'+self.__rtype+'"'+params_str+')'

    def __eq__(self,other):
        '''
        Check if two forward reaction rate coefficients are the same. They are same if the k values are equal.
        '''
        return self.k_forward() == other.k_forward()

    def __check_param_in(self, param_list):
        '''
        Returns dictionary of parameter values
        INPUT
        =====
        param_list: iterable, required
                    list of parameter namesof interest
        RETURNS
        =======
        Dictionary of parameters in provided list (returns None if mismatch present between parameter list and instance attributes)
        '''
        param_dict = dict()
        for param in param_list:
            if param in self.__params:
                param_dict[param] = self.__params[param]
            else:
                return None
        return param_dict

    def k_forward(self):
        '''
        Computes forward reaction rate coefficient

        NOTES
        ==========
        call corresponding private method to calculate k_values.

        RETURNS
        ==========
        k: (real number) Forward reaction rate coefficient
        '''
        if self.__rtype == "Constant":
            params = self.__check_param_in(['k'])
            if params != None:
                return self.__const(**params)
            else:
                raise ValueError('Insufficient Parameters')
        elif self.__rtype == "Arrhenius":
            params = self.__check_param_in(['A','E','T','R'])
            if params != None:
                return self.__arr(**params)
            else:
                raise ValueError('Insufficient Parameters')
        elif self.__rtype == "modifiedArrhenius":
            params = self.__check_param_in(['A','b','E','T','R'])
            if params != None:
                return self.__mod_arr(**params)
            else:
                raise ValueError('Insufficient Parameters')
        else:
            raise NotImplementedError("Type not supported yet!")

    def __const(self,k):
        '''
        Return constant coefficient
        INPUTS
        =======
        k: float, required
           constant reaction rate coefficient
        RETURNS
        ========
        k: float
           constant reaction rate coefficient
        '''

        if k < 0:
            raise ValueError("k cannot be negative!")
        return k

    def __arr(self,A,E,T,R):
        '''
        Return Arrhenius reaction rate coefficient
        INPUTS
        =======
        A: float, required
           Arrhenius prefactor  A. A  is strictly positive
        E: float, required
           Activation energy
        T: float, required
           Temperature. T must be positive (assuming a Kelvin scale)
        R: float, required
           Ideal gas constant
        RETURNS
        ========
        k: float
           Arrhenius reaction rate coefficient
        NOTES
        ========
        R = 8.314 is the default ideal gas constant.
        It should be positive (except to convert units)
        '''

        if A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

        return A * np.exp(-E / (R * T))

    def __mod_arr(self,A,b,E,T,R):
        '''
        Return modified Arrhenius reaction rate coefficient
        INPUTS
        =======
        A: float, required
           Arrhenius prefactor  A. A  is strictly positive
        b: float, required
           Modified Arrhenius parameter b. b must be real
        E: float, required
           Activation energy
        T: float, required
           Temperature. T must be positive (assuming a Kelvin scale)
        R: float, required
           The ideal gas constant
        RETURNS
        ========
        k: Modified Arrhenius reaction rate coefficient
        NOTES
        ========
        R=8.314 is the default ideal gas constant
        '''
        if A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

        return A * (T**b) * np.exp(-E / (R * T))
