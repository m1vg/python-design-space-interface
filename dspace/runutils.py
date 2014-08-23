'''Utility simplifying routine analysis by system design space.

'''
from __future__ import division
import dspace
import dspace.plotutils
import matplotlib as mt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.widgets import Slider
from math import *

class SliderCallback(object):
    
    def __init__(self, pvals, parameter, redraw):
        self._pvals = pvals
        self._name = parameter
        self._redraw = redraw
        
    def __call__(self, val):
        pvals = self._pvals
        pvals[self._name] = 10**val
        self._redraw()
            
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
                the name of an auxiliary variable. Default is no auxiliary 
                variables.
            
            constraints (list): A list of strings, where each string is
                a condition of the form 'lhs > 'rhs', or 'lhs < rhs'. The 
                constraints are added to the design space object.

            resolve_cycles (bool): A boolean value indicating if cycles should 
                be resolved automatcally. Default is False. 
                
            name (str): Used to identify a system.
                
            authors (list): A list of strings with the name of each author. 
                
            version (str): String representation of a version number. 
            
            get_parameters (int or str) : An int or str indicating the case
                number, subcase number or case signature marker for which an
                internal parameter set will be obtained.
            
            parameters (dict): A dictionary of the system parameters.
                         
            xaxis (str): Name of the variable on the x-axis of plots.
            
            yaxis (str): Name of the variable on the y-axis of plots.
            
            x_range (list): A list with the lower and upper bound of the 
                x-axis.
            
            y_range (list): A list with the lower and upper bound of the 
                y-axis.
            
            zlim (list): A list with the lower and upper bound of the 
                z-axis. This argument forces all heatmaps for concentrations,
                fluxes and arbitrary functions to have a fixed scale.
                
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
        constraints = None
        if 'auxiliary_variables' in options:
            aux = options['auxiliary_variables']
        if 'zlim' not in options:
            options['zlim'] = None
        eq = dspace.Equations(f_str, auxiliary_variables=aux)
        options.pop('equations')
        ds = dspace.DesignSpace(eq, **options)
        setattr(self, '_ds', ds)
        self._print_valid_cases(options)
        self._process_state(options)
        self._interactive_plot(options)
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
        if 'get_parameters' in options:
            case_id = options['get_parameters']
            case = self._ds(case_id)
            parameters = case.valid_interior_parameter_set(distance=1e3)
            options['parameters'] = parameters
            for i in pvals:
                print str(i) + ': ' + str(parameters[i])
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
        parameter_test = dict(options['parameters'])
        parameter_test[self._xaxis] = self._xrange
        parameter_test[self._yaxis] = self._yrange
        valid_test = self._ds.valid_cases(p_bounds=parameter_test)
        if len(valid_test) == 0:
            raise ValueError, 'Slice of parameter space does not produce valid cases' 

    def _print_valid_cases(self, options):
        if 'print_valid_cases' not in options:
            return
        if options['print_valid_cases'] is not True:
            return
        cases = self._ds.valid_cases()
        case_string = 'Valid Cases:\n'
        for i in cases:
            case = self._ds(i)
            case_string += str(i) + ': ' + case.signature + '\n'
        print case_string
    
    def _redraw_interactive(self):
        
        ax = self._plot_ax
        c_ax = self._plot_c_ax
        cdict = self._cdict
        ax.clear()
        c_ax.clear()
        ## print self._pvals['a1']
        cdict.update(self._ds.draw_2D_slice(ax, self._pvals,
                                            self._xaxis, self._yaxis,
                                            self._xrange, self._yrange,
                                            included_cases=self._included_cases,
                                            colorbar=False,
                                            color_dict=cdict))
        self._ds.draw_region_colorbar(c_ax, cdict)
        
    def _interactive_plot(self, options):
        
        if 'plot_interactive' not in options:
            return
        slider_dict = options['plot_interactive']
        if isinstance(slider_dict, dict) is False:
            raise TypeError, 'Interactive plot requires a dictionary of parameter : range pairs.'
        previous = plt.isinteractive()
        plt.ioff()
        number_of_sliders = len(slider_dict)
        slider_block = 0.03*number_of_sliders
        fig = plt.figure()
        plt.clf()
        self._plot_ax = plt.axes([0.1, 0.2+slider_block, 0.6, 0.7-slider_block])
        self._plot_c_ax = plt.axes([0.8, 0.2+slider_block, 0.1, 0.7-slider_block])
        self._cdict = dict()
        j = 0
        sliders = list()
        callbacks = list()
        for i in slider_dict:
            on_change = SliderCallback(self._pvals, i, self._redraw_interactive)
            slider_ax = plt.axes([0.1, 0.1+j*0.03, 0.8, 0.02])
            slider = Slider(slider_ax, i, 
                            log10(slider_dict[i][0]), log10(slider_dict[i][1]), 
                            valinit=log10(self._pvals[i]), color='#AAAAAA')
            j += 1
            sliders.append(slider)
            callbacks.append(on_change)
        
        for j in xrange(len(sliders)):
            sliders[j].on_changed(callbacks[j])
        
        plt.show()
        plt.interactive(previous)
        
    def _plot_designspace(self, options):
        
        if 'plot_designspace' not in options:
            return
        if options['plot_designspace'] is not True:
            return
        ## if ds.number_of_cases > 1e5 and self._included_cases is not None:
        ##     colors = dict()
        ##     j = 0
        ##     for i in self._included_cases:
        ##         case = self._ds(i)
        ##         colors[i] = cm.hsv(j/len(self._included_cases))
        ##         case.draw_2D_slice(plt.gca(), self._pvals,
        ##                            self._xaxis, self._yaxis,
        ##                            self._xrange, self._yrange,
        ##                            fc=colors[i], ec='none'
        ##                            )
        ##         j += 1
        ##     c_ax,kw=mt.colorbar.make_axes(plt.gca())
        ##     c_ax.set_aspect(15)
        ##     self._ds.draw_region_colorbar(c_ax, colors)
        ##     return
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
        ## if ds.number_of_cases > 1e5 and self._included_cases is not None:
        ##     colors = dict()
        ##     j = 0
        ##     for i in self._included_cases:
        ##         case = self._ds(i)
        ##         colors[i] = cm.hsv(j/len(self._included_cases))
        ##         case.draw_2D_slice(plt.gca(), self._pvals,
        ##                            self._xaxis, self._yaxis,
        ##                            self._xrange, self._yrange,
        ##                            fc=colors[i], ec='none'
        ##                            )
        ##         j += 1
        ##     c_ax,kw=mt.colorbar.make_axes(plt.gca())
        ##     c_ax.set_aspect(15)
        ##     self._ds.draw_region_colorbar(c_ax, colors)
        ##     return
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
                                         included_cases=self._included_cases,
                                         zlim=options['zlim'])
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
                                         included_cases=self._included_cases,
                                         zlim=options['zlim'])
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
                                         included_cases=self._included_cases,
                                         zlim=options['zlim'])
            ax.set_title(r'f:='+fn+' plot')
        
    
    
    
                              
                 
                 
