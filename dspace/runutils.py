'''Utility simplifying routine analysis by system design space.

'''

import dspace
import dspace.plotutils
import matplotlib.pyplot as plt

class Input(object):
    
    def __init__(self, 
                 equations,
                 **kwargs):
        ''' Initializes a new object with the input parameters for a routine
            analysis.
        
        The object incorporates the  input parameters to construct a system
        design space and specify the analysis that is performed.  The analyses
        that are performed involve routine functions, such as plotting a 2D
        design space,  steady state concentrations, steady state fluxes,
        stability information, etc.  The input arguments are a set
        of keyword arguments, that affect the output behavior.
                 
        Args:
            equations: A list of equations in string format defining the 
                system to analyze.
        
        Kwargs:
            auxiliary_variables (list): A list of strings, where each string is
                the name of an auxiliary variable. Default is noe auxiliary 
                variables.
            
            resolve_cycles (bool): A boolean value indicating if cycles should 
                be resolved automatcally. Default is False. 
                
            name (str): Used to identify a system.
                
            authors (list): A list of strings with the name of each author. 
                
            version (str): String representation of a version number. 
            
            parameters (dict): A dictionary of the system parameters.
                         
            xaxis (str): Name of the variable on the x-axis of plots.
            
            yaxis (str): Name of the variable on the y-axis of plots.
            
            x_range (list): A list with the lower and upper bound of the 
                x-axis.
            
            y_range (list): A list with the lower and upper bound of the 
                y-axis.
                
            centered_axes (bool): A boolean value indicating if axes ranges 
                should be centered around the x-axis and y-axis values
                specified in the parameter set. Default is False.
            
            print_valid_cases (bool): Specifies if the systems valid cases 
                will be printed. Default is False.

            print_parameters (bool): Specifies if the systems parameters 
                will be printed. Default is False.
            
            plot_designspace (bool): Specifies if the system design space 
                will be plotted. Default is False.
            
            plot_steady_states (list): A list of concentration variables to
                draw. Default is an empty set of plots to draw.
                
            plot_fluxes (list): A list of state variables for which their 
                steady state flux will be plotted. Default is an empty set of 
                plots to draw.
            
            plot_stability (bool): Specifies if the local stability will be 
                plotted. Default is False.
                
            draw_cases (list): A list of cases that will be the only ones
                drawn.  If all cases should be drawn, value should be None.
            
        '''
        options = {}
        options.update(kwargs)
        options.update({'equations':equations})
        self._main(options)
        
    def _main(self, options):
        f_str = options['equations']
        aux = []
        if 'auxiliary_variables' in options:
            aux = options['auxiliary_variables']
        eq = dspace.Equations(f_str, auxiliary_variables=aux)
        options.pop('equations')
        ds = dspace.DesignSpace(eq, **options)
        setattr(self, '_ds', ds)
        self._print_valid_cases(options)
        self._process_state(options)
        self._plot_designspace(options)
        self._plot_steady_states(options)
        self._plot_fluxes(options)
        self._plot_stability(options)
        self._plot_functions(options)
        plt.show()
    
    def _process_state(self, options):
        pvals = dspace.VariablePool(names=self._ds.independent_variables)
        xaxis = None
        yaxis = None
        x_range = None
        y_range = None
        setattr(self, '_pvals', pvals)
        setattr(self, '_xaxis', xaxis)
        setattr(self, '_yaxis', yaxis)
        setattr(self, '_xrange', x_range)
        setattr(self, '_yrange', y_range)
        setattr(self, '_included_cases', None)
        try:
            pvals.update(options['parameters'])
        except:
            pass
        centered_axes = False
        try:
            centered_axes = options['centered_axes']
        except:
            pass
        try:
            self._xaxis = options['xaxis']
            self._yaxis = options['yaxis']
            x_range = options['x_range']
            y_range = options['y_range']
            if centered_axes is False:
                self._xrange = (min(x_range), max(x_range))
                self._yrange = (min(y_range), max(y_range))
            else:
                self._xrange = (pvals[self._xaxis]*min(x_range), 
                                pvals[self._xaxis]*max(x_range))
                self._yrange = (pvals[self._yaxis]*min(y_range), 
                                pvals[self._yaxis]*max(y_range))
        except:
            pass
        try:
            self._included_cases=options['draw_cases']
        except:
            pass

    
    def _print_valid_cases(self, options):
        if 'print_valid_cases' not in options:
            return
        if options['print_valid_cases'] is not True:
            return
        cases = self._ds.valid_cases()
        case_string = 'Valid Cases:\n'
        for i in self._ds(cases):
            case_string += str(i.case_number) + ': ' + i.signature + '\n'
        print case_string

    def _plot_designspace(self, options):
        
        if 'plot_designspace' not in options:
            return
        if options['plot_designspace'] is not True:
            return
        fig = plt.figure()
        plt.clf()
        ax = plt.gca()
        self._ds.draw_2D_slice(plt.gca(), self._pvals,
                               self._xaxis, self._yaxis,
                               self._xrange, self._yrange,
                               included_cases=self._included_cases)
        ax.set_title('Design space plot')
        
    def _plot_stability(self, options):
        
        if 'plot_stability' not in options:
            return
        if options['plot_stability'] is not True:
            return
        fig = plt.figure()
        plt.clf()
        ax = plt.gca()
        self._ds.draw_2D_positive_roots(plt.gca(), self._pvals,
                                        self._xaxis, self._yaxis,
                                        self._xrange, self._yrange,
                                        included_cases=self._included_cases)
        ax.set_title('Stability plot')
        
    def _plot_steady_states(self, options):
        
        if 'plot_steady_states' not in options:
            return
        resolution=100
        try:
            resolution = options['resolution']
        except:
            pass
        for dependent in options['plot_steady_states']: 
            if dependent not in self._ds.dependent_variables:
                raise NameError, 'No variable named: ' + dependent
            fig = plt.figure()
            plt.clf()
            ax = plt.gca()
            self._ds.draw_2D_ss_function(plt.gca(), 'log('+dependent +')',
                                         self._pvals,
                                         self._xaxis, self._yaxis,
                                         self._xrange, self._yrange,
                                         resolution=resolution,
                                         log_linear=False,
                                         included_cases=self._included_cases)
            ax.set_title('[$'+dependent+'$] plot')
        
    def _plot_fluxes(self, options):
        
        if 'plot_fluxes' not in options:
            return
        resolution=100
        try:
            resolution = options['resolution']
        except:
            pass
        for dependent in options['plot_fluxes']: 
            if dependent not in self._ds.dependent_variables:
                raise NameError, 'No variable named: ' + dependent
            fig = plt.figure()
            plt.clf()
            ax = plt.gca()
            self._ds.draw_2D_ss_function(plt.gca(), 'log(V_'+dependent +')',
                                         self._pvals,
                                         self._xaxis, self._yaxis,
                                         self._xrange, self._yrange,
                                         resolution=resolution,
                                         log_linear=False,
                                         included_cases=self._included_cases)
            ax.set_title(r'$V_{'+dependent+'}$ plot')
            
    def _plot_functions(self, options):
        
        if 'plot_functions' not in options:
            return
        resolution=100
        log_linear=False
        try:
            resolution = options['resolution']
        except:
            pass
        try:
            log_linear = options['log_linear']
        except:
            pass
        for fn in options['plot_functions']: 
            fig = plt.figure()
            plt.clf()
            ax = plt.gca()
            self._ds.draw_2D_ss_function(plt.gca(), fn,
                                         self._pvals,
                                         self._xaxis, self._yaxis,
                                         self._xrange, self._yrange,
                                         resolution=resolution,
                                         log_linear=log_linear,
                                         included_cases=self._included_cases)
            ax.set_title(r'f:='+fn+' plot')
        
    
    
    
                              
                 
                 
