import numpy as np
import pandas as pd
from .ChemSolver import ChemSolver
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import io
import os
import urllib, base64
import re

class ChemViz:
    '''
    The ChemViz module is for visualization of the ChemSolver result.
    '''
    def __init__(self, chemsol):
        '''
        INPUT
        =====
        chem: chemkin object, required
        '''
        self.chemsol = chemsol
        if self.chemsol._sol!=True:
            raise ValueError('ChemSolver object must be solved before passed in ChemViz!')

    def plot_network(self):
        '''
        returns a figure, a place holder implementation
        '''
        plt.ioff()
        fig = plt.figure()
        plt.plot(range(10, 20))
        return fig

    def plot_time_series(self):
        '''
        returns a series
        '''
        plt.ioff()
        fig = plt.figure()
        plt.plot(range(10, 20))
        return fig
    
    def __species_encoding(self,species):
        return [re.sub(r"(\d)",r"<sub>\1</sub>",x) for x in species]
    
    def __conc_encoding(self,conc_array,species_array):
        return ''.join(["<p>"+str(y)+":"+" "+str(x)+"</p>" for x,y in zip(conc_array,species_array)])
        
    def html_report(self,file_name):
        # check file name
        if not ('.html' == file_name[-5:]):
            raise ValueError('The filename suffix must be .html.')
            
        # get base64 encoded string
        buf = io.BytesIO()
        self.plot_network().savefig(buf,format='png')
        buf.seek(0)
        network_str = urllib.parse.quote(base64.b64encode(buf.read()))
        buf = io.BytesIO()
        self.plot_time_series().savefig(buf,format='png')
        buf.seek(0)
        time_series_str = urllib.parse.quote(base64.b64encode(buf.read()))

        #read html template
        here = os.path.abspath(os.path.dirname(__file__))
        template_f = os.path.join(here, 'data/template.html')
        if not os.path.isfile(template_f):
            raise ValueError('Template does not exist. Please reinstall')
        with open(template_f,'r') as f:
            template_str = f.read()

        # add plot
        template_str = template_str.replace("$base64_network$",network_str)
        template_str = template_str.replace("$base64_timeseries$",time_series_str)
        # add temperature, time
        template_str = template_str.replace("$temperature$",str(self.chemsol.T))
        template_str = template_str.replace("$end_time$",str(self.chemsol._t[-1]))
        
        # add concentration
        species = self.__species_encoding(self.chemsol.chem.species)
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
            for specie in reaction.split():
                if specie.find('=]')!=-1:
                    if see_equal_sign:
                        raise ValueError("Equation format is incorrect: "+str(reaction))
                    see_equal_sign = True
                    if self.chemsol.chem.reversible[reaction_cnt]:
                        cur_reaction_str += ' ⇌'
                    else:
                        cur_reaction_str += ' ='
                elif specie=='+':
                    cur_reaction_str += ' +'
                else:
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