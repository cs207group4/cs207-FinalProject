import numpy as np
import pandas as pd
from .ChemSolver import ChemSolver
import matplotlib.pyplot as plt
import io
import os
import urllib, base64

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
        plt.plot(range(10, 20))
        fig = plt.gcf()
        return fig

    def plot_time_series(self):
        '''
        returns a series
        '''
        plt.plot(range(10, 20))
        fig = plt.gcf()
        return fig

    def html_report(self,file_name):

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
        template_str.replace("$base64_network$",network_str)
        template_str.replace("$base64_timeseries$",time_series_str)
        if not ('.html' == file_name[-5:]):
            raise ValueError('The filename suffix must be .html.')
        with open(file_name,'w') as f:
            f.write(template_str)