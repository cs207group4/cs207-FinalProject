class ReactionCoeffs:
    """
    A class for computating reaction coefficients

    Initialize by calling reaction_coeffs(type, {params})
    Set params by calling set_params({params})
    Calculate k value by calling kval()

    EXAMPLES
    =========
    >>> rc = ReactionCoeffs('Constant', k = 1e3)
    >>> rc.kval()
    1000.0
    >>> rc = ReactionCoeffs('Arrhenius', A = 1e7, E=1e3)
    >>> rc.set_params(T=1e2)
    >>> rc.kval()
    3003549.0889639612
    """

    def __init__(self, type, **kwargs):
        """
        INPUT
        =====
        type: string, required
              Type of the reaction rate coefficient of interest (e.g. "Constant", "Arrhenius", etc)
        """
        if type not in ["Constant","Arrhenius","modifiedArrhenius"]:
            raise NotImplementedError('Rate coefficient type not supported yet! Must be constant, Arrhenius, or modified Arrhenius.')
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
        """
        Check if two reaction rate coefficients are the same. They are same if the k values are equal.
        """
        return self.kval() == other.kval()

    def __check_param_in(self, param_list):
        """
        Returns dictionary of parameter values

        INPUT
        =====
        param_list: iterable, required
                    list of parameter namesof interest

        RETURNS
        =======
        Dictionary of parameters in provided list (returns None if mismatch present between parameter list and instance attributes)
        """
        param_dict = dict()
        for param in param_list:
            if param in self.__params:
                param_dict[param] = self.__params[param]
            else:
                return None
        return param_dict

    def kval(self):
        """
        Computes reaction coefficient

        NOTES
        ==========
        call corresponding private method to calculate k_values.
        Easy to extend... Just write a new private method for new type of coeffs.

        RETURNS
        ==========
        k: (real number) Reaction coefficient
        """
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
        """
        Return constant coefficient

        INPUTS
        =======
        k: float, required
           constant reaction rate coefficient

        RETURNS
        ========
        k: float
           constant reaction rate coefficient

        """

        if k < 0:
            raise ValueError("k cannot be negative!")
        return k

    def __arr(self,A,E,T,R):
        """
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
        R = 8.314 is the ideal gas constant. It should never be changed (except to convert units)

        """
        # R should never be changed
        #for the first milestone we assume consistent units throughout
        if R!=8.314:
            raise ValueError("R should not be changed.")
        if A <= 0:
            raise ValueError("A must be positive!")
        if T <= 0:
            raise ValueError("T must be positive!")

        # On overflow, return np.inf
        try:
            return A*np.exp(-E/(R*T))
        except OverflowError:
            print("WARNING: Overflow error")
            return np.inf

    def __mod_arr(self,A,b,E,T,R):
        """
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
        R=8.314 is the ideal gas constant. It should never be changed (except to convert units)
        """
        # R should never be changed
        if A <= 0:
            raise ValueError("A must be positive!")
        if T <= 0:
            raise ValueError("T must be positive!")
        if not isinstance(b, (float, int)):
            raise ValueError("b must be real!")

        # On overflow, return np.inf
        try:
            return A*(T**b)*np.exp(-E/(R*T))
        except OverflowError:
            print("WARNING: Overflow error")
            return np.inf
