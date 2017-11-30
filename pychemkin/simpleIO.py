import pickle

class simpleIO:
    '''
    The simpleIO module wraps around pickle's IO methods to provide a simple way to save and load 
    python objects. Users can load the saved pychemkin objects (such as ChemSolver or ChemViz) without
    explicitly preloading dependent objects (such as chemkin)
    
    METHODS
    ========
    After specifying file_name (including path), user could call:
     - to_pickle(obj): method to save a python object to python-specific pickle format
     - read_pickle(): method to load a python object from a pickle file
    '''
    def __init__(self, file_name):
        '''
        INPUT
        =====
        file_name: string, required
            The name (including path) which user wants to save a python object to or load a python
            object from. The file extension must be '.pkl'; otherwise, a ValueError would be raised
        '''
        file_name = file_name.strip()
        if file_name[-4:] != '.pkl':
            raise ValueError('Only ".pkl" file is supported.')
        self.file_name = file_name
        
    def to_pickle(self, obj):
        '''
        INPUT
        =====
        obj: python object, required
        '''
        with open(self.file_name, 'wb') as output:
            pickle.dump(obj, output, protocol=pickle.HIGHEST_PROTOCOL)
        return self
    
    def read_pickle(self):
        with open(self.file_name, 'rb') as input_:
            obj = pickle.load(input_)
        return obj