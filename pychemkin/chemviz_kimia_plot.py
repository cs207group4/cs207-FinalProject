class ChemViz: 
    
    """
    This will allow users to visualize species’ reaction rates or concentrations over time.
    
    Input: 
    - file is a string or a list of strings corresponding the the csv or HDF5 file(s) 
        containing the simulation data to be plotted ChemSolver.solve.
    - species is a list of species to be plotted
    - LowT specifies the start time for the time series to be plotted 
        (which can be a subset of the simulation run).
    - HighT specifies the end time for the time series to be plotted 
        (which can be a subset of the simulation run). 
    - yaxis is the quantity to be plotted: reaction rate or concentration

    """
        
    def __init__(self):
        
        self.file = None
        self.species = None
        ## species will be a list
        self.T = None
        ## will also be a list/range
        self.simulation_file = None
        self.yaxis= None
        
    
    def readinfile(self, file):

        filename = str(file)
        last_3_letters = filename.strip()[-3:]
        if last_3_letters == 'csv': 
            self.df = pd.read_csv(str(file))
        elif last_3_letters == 'df5': 
            self.df =  pd.HDFStore(str(file))
        elif last_3_letters == '.h5': 
            self.df = pd.HDFStore(str(file))
        else: 
            raise ValueError("file must be CSV or H5/HDF5")
    
        
    def filterByTemperature(self, LowT, HighT):
        return self.df[(self.df['t']<=HighT) & (self.df['t']>=LowT)]

        
    def plot_time_series(self, file, species, LowT, HighT, yaxis, outputfile=None):
        
        #read in file
        self.readinfile(str(file))
        
        ## Concentration or reaction rate
        if type(yaxis) != str:
            raise TypeError("yaxis must be string")
        ## Concentration 
        lowyaxis = yaxis.lower()
        if lowyaxis == "concentration": 
            concentration_cols = [col for col in self.df.columns if 'Concentration' in col] + ['t']
            self.df = self.df[concentration_cols]
        ## reaction rate
        lowyaxis = yaxis.lower()
        if lowyaxis == "reactionrate": 
            rxn_cols = [col for col in self.df.columns if 'Reaction_rate' in col] + ['t']
            self.df = self.df[rxn_cols]
        
    
        #specify for what temp
        self.df = self.filterByTemperature(LowT, HighT)
        
        #rename the columns of the df
        columns = (self.df.columns)
        for i,v in enumerate(columns):
            
            column_new_name = (columns[i].split('-')[0])
            self.df = self.df.rename(columns={str(v): str(column_new_name)})
        print(self.df.columns)
            
        #Specify which species
        listofsp = list(species)
        self.dfofsp = self.df[listofsp]

        
        #plot
        ### what to do with this below
        %matplotlib inline
        plt.figure(figsize = (10,5))
        ###?????? NEED TO WORK ON COLORZ
        list_of_colors = ((33, 55, 66,))
        
        title = "plot of " + str(species)

        plt.title(str((title)))  
        
        for i,v in enumerate(species):
            plt.plot(self.df['t'],self.dfofsp[str(v)])
            
        plt.legend()
