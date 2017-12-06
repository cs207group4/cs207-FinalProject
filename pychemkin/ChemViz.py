import io
import os
import urllib, base64
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from .chemkin import chemkin
from .ChemSolver import ChemSolver

class ChemViz:
    '''
    The ChemViz module is for visualization of ODE solutions obtained by the ChemSolver module

    INSTANCE METHODS and ATTRIBUTES
    ========
     - plot_network(yaxis, species = None, tmin = 0, tmax = None, outputfile=None): plots reaction diagram
     - plot_time_series(yaxis, species = None, tmin = 0, tmax = None, outputfile=None): plots a single time series of concentration or reaction rates
     - plot_gridtime_series(yaxis, species = None, tmin = 0, tmax = None, outputfile=None): plots a grid of time series for temperature/initial concentration combos
     - html_report(file_name): generate an HTML report and save it to file_name

    See method doc for additional details

    EXAMPLE
    =========
    >>> chem = chemkin("tests/test_xml/rxns.xml")
    Finished reading xml input file
    >>> y0 = np.ones(len(chem.species))
    >>> T = 300 #temperature in K
    >>> t1 = 0.05 #Last timepoint for which to solve ODE
    >>> dt = 2e-4 #time step for solutions to be returned
    >>> cs = ChemSolver(chem).solve(y0, T, t1, dt, method='lsoda')
    >>> ChemViz(cs).html_report('report.html')
    '''
    def __init__(self, chemsol):
        '''
        INPUT
        =====
        chemsol: solved ChemSolver object, required
        '''
        self.chemsol = chemsol
        if self.chemsol._sol!=True and self.chemsol.grid_result is None:
            raise ValueError('ChemSolver object passed to ChemViz needs to have the solutions attribute set or grid_result attribute set.')
        if self.chemsol._sol:
            self.df = chemsol.to_df()
            self.end_time = self.df['t'][-1:].values[0]

    def _filterByTime(self, df, tmin, tmax):
        """
        Select ODE solutions in a particular time interval

        INPUT
        =====
        df: Pandas dataframe, required
            Dataframe holding time series data
        tmin: Float, required
              Start of desired interval
        tmax: Float, required
              End of desired interval
        """
        end_time = df['t'][-1:].values[0]
        if tmin < 0:
            raise ValueError('Starting time has to be >=0 s')
        if tmax > end_time: #last timepoint for which data is available
            print("WARNING: Last time point is {} s".format(end_time))
        return df[(df['t']<=tmax) & (df['t']>=tmin)]


    def _subplot_network(self, fig, nrow, ncol, panel, timepoint):
        """
        Plots reaction diagram for a single timepoint

        INPUT
        =====
        fig: Matplotlib Figure instance, required
        nrow: Integer, required
              Number of rows desired for network plot
        ncol: Integer, required
              Number of columns desired for network plot
        panel: Integer, required
              Panel number of subplot
        timepoint: Float, required
              Timepoint for which subplot is being made

        RETURNS
        =======
        Matplotlib Ax instance with the network subplot information
        """
        if timepoint < 0 or timepoint > self.end_time:
            raise ValueError('Time point must be between 0 and {} seconds'.format(self.end_time))

        T = self.chemsol.T #temperature in K

        chemkinobj = self.chemsol.chem
        species = chemkinobj.species
        N, M = chemkinobj.nu_react.shape #N species and M reactions
        #interpolate concentrations
        interp_conc = interp1d(self.df['t'].values, self.df.iloc[:, 1:1+N].values, axis = 0)
        #interpolate reaction rates
        interp_rxnrate = interp1d(self.df['t'].values, self.df.iloc[:, 1+N:].values, axis = 0)

        conc = interp_conc(timepoint)
        rxnrate = interp_rxnrate(timepoint)
        _,kf, kb = chemkinobj._progress_rate_init(np.ones(N), T) #putting in a dummy array for concentration just to get reaction rate coefficients
        #rescale the reaction coefficients in order to set connecting line weights
        mincoeff = np.min(np.concatenate((kf,kb)))
        normed_kf = np.log10(kf/mincoeff)
        if len(kb)>0:
            normed_kb = np.log10(kb/mincoeff)
        else:
            normed_kb = None


        angle = np.arange(0, 2*np.pi,2*np.pi/N)
        ax = fig.add_subplot(nrow, ncol, panel,aspect = 'equal')

        normalized_conc = conc/np.max(conc)
        specie_position = []


        for i, specie in enumerate(species):
            specie_position.append((30*np.cos(angle[i]),30*np.sin(angle[i])))

        fwd_rxn = []
        back_rxn = []
        for i in range(M):
            fwd_rxn.append(np.where(chemkinobj.nu_react[:,i]>0)[0])
            back_rxn.append(np.where(chemkinobj.nu_prod[:,i]>0)[0])
        reversible_indices = np.where(chemkinobj.reversible)[0]
        #plot lines connecting species that react with one another
        for i in range(M):
            if chemkinobj.reversible[i]:
                if len(back_rxn[i])==2:
                    specie1, specie2 = back_rxn[i]
                    x1,y1 = specie_position[specie1]
                    x2,y2 = specie_position[specie2]
                    ax.arrow(x1,y1,
                        x2-x1,y2-y1,
                        fc="b",
                        ec="b",
                        width = 0.1 + 0.1*normed_kb[np.where(reversible_indices==i)[0][0]],
                        head_width=0,
                        head_length=0,
                        alpha = 0.5)
                else:
                    #account for cases where specie reacts with itself
                    specie1 = specie2 = back_rxn[i][0]
                    loop = Circle(xy = np.array(specie_position[specie1])*1.15,
                                  radius = 3.5,
                                  fill = False,
                                  color = 'b')
                    ax.add_artist(loop)

            if len(fwd_rxn[i])==2:
                specie1, specie2 = fwd_rxn[i]
                x1,y1 = specie_position[specie1]
                x2,y2 = specie_position[specie2]
                ax.arrow(x1,y1,
                    x2-x1,y2-y1,
                    fc="b",
                    ec="b",
                    width = 0.1+0.1*normed_kf[i],
                    head_width=0,
                    head_length=0,
                    alpha = 0.5)

            else:
                specie1 = specie2 = fwd_rxn[i][0]
                loop = Circle(xy = np.array(specie_position[specie1])*1.15,
                                  radius = 3.5,
                                  fill = False,
                                  color = 'b')
                ax.add_artist(loop)
        #plot bubbles with species names
        for i, specie in enumerate(species):
            if rxnrate[i]<0:
                color = '#b30000'
            else:
                color = '#007a99'
            outer_circ= Circle(xy = specie_position[i],radius = 3+4*normalized_conc[i], fill = True, color = color)
            ax.add_artist(outer_circ)
            inner_circ= Circle(xy = specie_position[i],radius = 3, color = 'white')
            ax.add_artist(inner_circ)
            ax.text(30*np.cos(angle[i]),
                    30*np.sin(angle[i]),
                    self.__make_tex_string(specie),
                    color = 'black',
                    size = 14,
                    horizontalalignment = 'center',
                    verticalalignment = 'center')



        ax.text(0, 40, "$t =$ %.3e s" % timepoint, horizontalalignment = 'center', size = 30)
        ax.axis('off')
        ax.set_xlim(xmin = -45, xmax = 45)
        ax.set_ylim(ymin = -45, ymax = 45)

        return ax

    def plot_network(self, timepoints, outputfile = None, figsize = (7,14)):
        """
        Make diagram of reaction network at one or more timepoints

        INPUT
        =====
        timepoints: Array of floats, required
                    Time points for which a reaction diagram is of interest
        outfile: String, required
                Name of png file to save to, if desired
        """
        if self.chemsol._sol!=True:
            raise ValueError('Check that ChemSolver object has solutions attribute set')

        plt.ioff()
        figsize = figsize
        fig = plt.figure(figsize = figsize)
        for i, t in enumerate(timepoints):
            self._subplot_network(fig, len(timepoints), 1, i+1, t)
        plt.tight_layout()

        if outputfile:
            plt.savefig(outputfile)
        return fig

    def _plot_df(self, dataframe, fig, nrow, ncol, panel, yaxis, T, species = None, tmin = 0, tmax = None):
        """
        Helper function for making a plot of either species concentrations or reaction rates as a function of time for a single simulation from ChemSolver

        INPUT
        ====
        dataframe: pandas dataframe, required
               dataframe with column information about timesteps, concentration rates, and reaction rates
        fig: Matplotlib figure instance, required
        nrow: Integer, required
              Number of rows desired for network plot
        ncol: Integer, required
              Number of columns desired for network plot
        panel: Integer, required
              Panel number of subplot
        yaxis: String, required
               Parameter to be plotted ('reactionrate' or 'concentration')
        T: float, required
               Temperature at which system of reactions takes place
        species: List of strings, required
               Species to be plotted. Default is all species
        tmin: float, required
              First time point in plot. Default is 0s
        tmax: float, required
              Last time point in plot. Default is last time point in simulation
        """
        #make default species list all species
        if species is None:
            species = self.chemsol.chem.species
        #set default for time range to be entire time range solved for by ODE solver
        if tmax is None:
            tmax = self.end_time


        ## Concentration or reaction rate
        if type(yaxis) != str:
            raise TypeError("yaxis must be string")
        filtered_df = dataframe.copy()
        ## Concentration
        lowyaxis = yaxis.lower()
        if lowyaxis == "concentration":
            concentration_cols = [col for col in dataframe.columns if 'Concentration' in col] + ['t']
            filtered_df = filtered_df[concentration_cols]
        ## reaction rate
        elif lowyaxis == "reactionrate":
            rxn_cols = [col for col in dataframe.columns if 'Reaction_rate' in col] + ['t']
            filtered_df = filtered_df[rxn_cols]

        else:
            raise ValueError('values for yaxis must be either "concentration" or "reactionrate"')


        #specify time range
        filtered_df= self._filterByTime(filtered_df, tmin, tmax)

        #rename the columns of the df
        columns = (filtered_df.columns)
        for i,v in enumerate(columns):

            column_new_name = (columns[i].split('-')[0])
            filtered_df = filtered_df.rename(columns={v: column_new_name})

        #Specify which species
        listofsp = list(filtered_df.columns)
        filtered_df = filtered_df[listofsp]
        ax = fig.add_subplot(nrow, ncol, panel)
        title = "{} vs. time at T = {} K".format(lowyaxis, T)

        ax.set_title(title)

        for i,v in enumerate(species):

            ax.plot(filtered_df['t'],filtered_df[v], label = self.__make_tex_string(v), color = plt.cm.brg((i+1)/(len(species))))
        if panel==nrow*ncol:
            ax.legend(bbox_to_anchor=(1.25, 1.1))
        ax.set_xlim(xmin = tmin, xmax = tmax)
        ax.set_ylabel(lowyaxis)
        ax.set_xlabel('time (seconds)')
        return ax



    def plot_time_series(self, yaxis, species = None, tmin = 0, tmax = None, outputfile=None):
        """
        Make plot of either species concentrations or reaction rates as a function of time for a single simulation from ChemSolver

        INPUT
        ====
        yaxis: String, required
               Parameter to be plotted ('reactionrate' or 'concentration')
        species: List of strings, required
               Species to be plotted. Default is all species
        tmin: float, required
              First time point in plot. Default is 0s
        tmax: float, required
              Last time point in plot. Default is last time point in simulation
        outputfile: String, optional
             Name of png file to save to
        """
        if self.chemsol._sol!=True:
            raise ValueError('Check that ChemSolver object has solutions appropriately set.')

        plt.ioff()
        fig = plt.figure(figsize = (7,4))

        self._plot_df(self.df, fig, 1, 1, 1, yaxis, self.chemsol.T, species = species, tmin = tmin, tmax = tmax)
        plt.subplots_adjust(right = 0.8)
        if outputfile:
            plt.savefig(outputfile)

        return fig

    def plot_gridtime_series(self, yaxis, species = None, tmin = 0, tmax = None, outputfile=None):
        """
        Make plot of either species concentrations or reaction rates as a function of time for grid simulations from ChemSolver

        INPUT
        ====
        yaxis: String, required
               Parameter to be plotted ('reactionrate' or 'concentration')
        species: List of strings, required
               Species to be plotted. Default is all species
        tmin: float, required
              First time point in plot. Default is 0s
        tmax: float, required
              Last time point in plot. Default is last time point in simulation
        outputfile: String, optional
             Name of png file to save to
        """
        plt.ioff()
        if self.chemsol.grid_condition is None or self.chemsol.grid_result is None:
            raise ValueError('No grid solutions have been computed yet!')
        num_temps = len(self.chemsol.grid_condition[0])
        num_conc = len(self.chemsol.grid_condition[1])
        fig = plt.figure(figsize = (6*num_temps,6*num_conc))
        for i, Temp in enumerate(self.chemsol.grid_condition[0]):
            for j, conc in enumerate(self.chemsol.grid_condition[1]):
                df = self.chemsol._gridsol_to_df(Temp, self.chemsol.grid_result[Temp][j])
                self._plot_df(df, fig, num_temps, num_conc, i*num_conc+j+1, yaxis, Temp, species = species, tmin = tmin, tmax = tmax)
        plt.subplots_adjust(right = 0.75, hspace = 0.25, wspace = 0.25)
        if outputfile:
            plt.savefig(outputfile)

        return fig

    def __make_tex_string(self, specie):
        '''
        Helper function to format the species name with proper underscores in LaTeX

        INPUT
        =====
        specie: string, required
             Specie name (as printed in original xml input file)

        RETURNS
        =======
        LaTeX-formatted string for specie name

        '''
        numeral_list = ['0','1','2','3','4','5','6','7','8','9']
        newstring = ""
        for i, char in enumerate(specie):
            if char in numeral_list:
                if specie[i-1] not in numeral_list:
                    newstring+='$_{}'.format(char)
                else:
                    newstring+=char
                if i==len(specie)-1 or specie[i+1] not in numeral_list:
                    newstring+='$'
            else:
                newstring+=char
        return newstring


    def __species_encoding(self,species):
        '''Helper function to format the species names with proper underscores in HTML

        INPUT
        =====
        species: List of strings, required
                 List of species names
        RETURNS
        =======
        List of HTML strings representing the names of each species
        '''
        return [re.sub(r"(\d)",r"<sub>\1</sub>",x) for x in species]

    def __conc_encoding(self,conc_array,species_array):
        '''helper function to format the concentration string in HTML

        INPUT
        =====
        conc_array: Array of floats, required
                Array of species concentrations

        species_array: List of strings, required
                List of species names corresponding to the concentration vector
        RETURNS
        =======
        List of HTML strings representing the concentration of each species
        '''
        if len(conc_array)!=len(species_array):
            raise ValueError('conc_array and species_array must match in size')
        return ''.join(["<p>"+str(y)+":"+" "+str(x)+"</p>" for x,y in zip(conc_array,species_array)])

    def html_report(self,file_name):
        '''Generate an HTML report of the reaction system dynamics.
        It contains:
            - the reactions system, conditions and coefficients
            - initial and ending concentrations
            - whether equilibrium has been reached
            - reaction network plot
            - reaction time series

        INPUT
        =======
        file_name: a string indicating the path where the html file should goto. Must end with .html

        RETURN
        =======
        No return. Function will save an HTML file with the file_name as the file name.

        '''

        # check file name
        if not ('.html' == file_name[-5:]):
            raise ValueError('The filename suffix must be .html.')

        # get base64 encoded string
        buf = io.BytesIO()
        self.plot_network(np.array([0, self.df['t'][-1:].values[0]])).savefig(buf,format='png')
        buf.seek(0)
        network_str = urllib.parse.quote(base64.b64encode(buf.read()))
        buf = io.BytesIO()
        self.plot_time_series('concentration').savefig(buf,format='png')
        buf.seek(0)
        concentration_str = urllib.parse.quote(base64.b64encode(buf.read()))
        buf = io.BytesIO()
        self.plot_time_series('reactionrate').savefig(buf,format='png')
        buf.seek(0)
        rxnrate_str = urllib.parse.quote(base64.b64encode(buf.read()))

        #read html template
        here = os.path.abspath(os.path.dirname(__file__))
        template_f = os.path.join(here, 'data/template.html')
        if not os.path.isfile(template_f):
            raise ValueError('HTMLTemplate does not exist. Please reinstall pychemkin.')
        with open(template_f,'r') as f:
            template_str = f.read()

        # add plot
        template_str = template_str.replace("$base64_network$",network_str)
        template_str = template_str.replace("$base64_concentration$",concentration_str)
        template_str = template_str.replace("$base64_rxnrate$",rxnrate_str)

        # add temperature, time
        template_str = template_str.replace("$temperature$",str(self.chemsol.T))
        template_str = template_str.replace("$end_time$",str(self.chemsol._t[-1]))

        # add concentration
        species = self.__species_encoding(self.chemsol.chem.species)
        conc_initial = self.chemsol._y[:,0]
        conc_final = self.chemsol._y[:,-1]
        conc_string = ''
        for i, specie in enumerate(self.chemsol.chem.species):
            conc_string += '<tr><td>{}</td><td>{}</td><td>{}</td></tr>'.format(species[i], conc_initial[i], conc_final[i])
        template_str = template_str.replace("$conc_string$",conc_string)


        template_str = template_str.replace("$ini_conc$",self.__conc_encoding(self.chemsol._y[:,0],species))
        template_str = template_str.replace("$end_conc$",self.__conc_encoding(self.chemsol._y[:,-1],species))
        template_str = template_str.replace("$is_equilibrium$",
                                           "<p style='color:blue'>System has reached equilibrium</p>"
                                            if self.chemsol.is_equilibrium()
                                           else "<p style='color:red'>System hasn't reached equilibrium yet</p>")

        # add reaction system
        reaction_str = ''
        reaction_cnt = 0
        reverse_reaction_cnt = 0

        for reaction in self.chemsol.chem.equations:
            reaction_str += '<tr>'
            cur_reaction_str = ''
            see_equal_sign = False
            specie_set = set()
            for specie in reaction.split():
                if specie.find('=]')!=-1:
                    if see_equal_sign:
                        raise ValueError("Equation format is incorrect: "+str(reaction))
                    see_equal_sign = True
                    if self.chemsol.chem.reversible[reaction_cnt]:
                        cur_reaction_str += ' &#8644;'
                    else:
                        cur_reaction_str += ' &#8594;'
                elif specie=='+':
                    cur_reaction_str += ' +'
                else:
                    if specie in specie_set:
                        # remove the previous +
                        cur_reaction_str = cur_reaction_str[:-2]
                        continue
                    specie_set.add(specie)
                    nu_react = self.chemsol.chem.nu_react[self.chemsol.chem.species.index(specie),reaction_cnt]
                    nu_prod = self.chemsol.chem.nu_prod[self.chemsol.chem.species.index(specie),reaction_cnt]
                    nu = nu_prod if see_equal_sign else nu_react
                    nu = str(nu) if nu!=1 else ''
                    cur_reaction_str += ' '+nu+ self.__species_encoding([specie])[0]
            reaction_str += '<td>' + cur_reaction_str.strip()+'</td>'
            reaction_str += '<td>' + str(self.chemsol.kf[reaction_cnt])+'</td>'
            kb = 'N/A'
            if self.chemsol.chem.reversible[reaction_cnt]:
                kb = str(self.chemsol.kb[reverse_reaction_cnt])
                reverse_reaction_cnt += 1
            reaction_str += '<td>' + kb + '</td>'
            reaction_str += '</tr>'
            reaction_cnt += 1
        template_str = template_str.replace("$reaction_system$",reaction_str)
        #save
        with open(file_name,'w') as f:
            f.write(template_str)
