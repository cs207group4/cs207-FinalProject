class ReactionCoeffs:
    ''' A class for reaction_coeffs
        
        Initialize by calling reaction_coeffs(type, {params})
        Set params by calling set_params({params})
        Get calculated k value by calling kval()
        
        EXAMPLES
        =========
        >>> rc = ReactionCoeffs('Constant', k = 1e3)
        >>> rc.kval() 
        1000.0
        >>> rc = ReactionCoeffs('Arrhenius', A = 1e7, E=1e3)
        >>> rc.set_params(T=1e2)
        >>> rc.kval() 
        3003549.0889639612
    '''

    def __init__(self, type, **kwargs):
        if type not in ["Constant","Arrhenius","modifiedArrhenius"]:
            raise NotImplementedError('Type Not Supported Yet!')
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
        return self.kval() == other.kval()

    def __check_param_in(self, param_list):
        '''
        return None if not all params avaliable
        else return a dict containing all avaliable params
        '''
        param_dict = dict()
        for param in param_list:
            if param in self.__params:
                param_dict[param] = self.__params[param]
            else:
                return None
        return param_dict

    def kval(self):
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
        ''' return constant coefficient
        
        INPUTS
        =======
        k: constant coefficient
        
        RETURNS
        ========
        k constant coefficient
        
        
        '''
        
        if k < 0:
            raise ValueError("k can not be negative!")
        return k

    def __arr(self,A,E,T,R):
        ''' return Arrhenius reaction rate coefficient
        
        INPUTS
        =======
        A: Arrhenius prefactor  A. A  is strictly positive
        E: Activation energy
        T: Temperature. T must be positive (assuming a Kelvin scale)
        R: Ideal gas constant
        
        
        RETURNS
        ========
        k: Arrhenius reaction rate coefficient
        
        NOTES
        ========
        R=8.314 is the ideal gas constant. It should never be changed (except to convert units)
        
        '''
        # R should never be changed
        if A <= 0:
            raise ValueError("A must be positive!")
        if T <= 0:
            raise ValueError("T must be positive!")
        
        # On overflow, return np.inf
        try:
            return A*np.exp(-E/(R*T))
        except OverflowError:
            return np.inf

    def __mod_arr(self,A,b,E,T,R):
        ''' return modified Arrhenius reaction rate coefficient
        
        INPUTS
        =======
        A: Arrhenius prefactor  A. A  is strictly positive
        b: Modified Arrhenius parameter b. b must be real
        E: Activation energy
        T: Temperature. T must be positive (assuming a Kelvin scale)
        R: The ideal gas constant
        
        
        RETURNS
        ========
        k: Modified Arrhenius reaction rate coefficient
        
        NOTES
        ========
        R=8.314 is the ideal gas constant. It should never be changed (except to convert units)

        
        '''
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
            return np.inf