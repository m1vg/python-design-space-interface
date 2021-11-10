import dspace
import dspace.plotutils
import dspace.display
#from InteractiveInput import load_widget
import cPickle as pickle

import numpy as np

from distutils.version import LooseVersion, StrictVersion

import IPython

if StrictVersion(IPython.__version__) < StrictVersion('4.0.0'):
    from IPython.html.widgets import interact, interactive, fixed
    from IPython.html.widgets import HTMLWidget as HTML
    from IPython.html.widgets import TabWidget as Tab
    from IPython.html.widgets import CheckboxWidget as Checkbox
    from IPython.html.widgets import ButtonWidget as Button
    from IPython.html.widgets import ContainerWidget as Box
    from IPython.html.widgets import TextWidget as Text
    from IPython.html.widgets import TextareaWidget as Textarea
    from IPython.html.widgets import DropdownWidget as Dropdown
    from IPython.html.widgets import RadioButtonsWidget as RadioButtons
    from IPython.html.widgets import PopupWidget as Popup
    from IPython.html.widgets import LatexWidget as Latex
    from IPython.html.widgets import FloatTextWidget as FloatText
    from IPython.html.widgets import ImageWidget as Image
    from IPython.html.widgets import IntSliderWidget as Slider
    VBox = Box
    HBox = Box
    old_ipython = True
else:
    from ipywidgets import *
    from popup import Popup
    old_ipython = False
    
from IPython.display import clear_output, display

import matplotlib as mt
import matplotlib.pyplot as plt

import cStringIO
import base64
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from subprocess import call, Popen, PIPE
from dspace.graphs.designspace_graph import GraphGenerator

def sort_eigenvalues(a, b):
    if a.real > b.real:
        return -1
    elif b.real > a.real:
        return 1
    else:
        return 0
    
def eigenvalue_compare(eigenvalues, component='real', rank=1):
    eig = [(eigenvalues.real[i], eigenvalues.imag[i]) for i in range(len(eigenvalues))]
    eig = sorted(eig, key = lambda i: i[0])
    if component == 'real':
        value = [i[0] for i in eig]
    else:
        value = [i[1] for i in eig]
    rank = min(rank, len(eig))
    return value[-rank]

class MakePlot(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'plot_data', VBox())
        setattr(self, 'title', None)
        setattr(self, 'caption', None)
        setattr(self, 'is_1D', False)
    
    @property
    def widget_types(self):
        widget_types = ['Design Space (interactive)',
                        'Design Space',
                        'Steady State Concentration',
                        'Steady State Flux',
                        'Steady State Function',
                        'Stability',
                        'Eigenvalues'
                        ]
        return widget_types
        
    def create_plot_widget(self):
        controller = self.controller
        options = self.controller.options
        if controller.ds is None:
            return
        xaxis = controller.defaults('xaxis')
        if xaxis is None:
            xaxis = controller.ds.independent_variables[0]
        yaxis = controller.defaults('yaxis')
        if yaxis is None:
            yaxis = controller.ds.independent_variables[1]
        center = controller.defaults('center_axes')
        range_x = controller.defaults('range_x')
        range_y = controller.defaults('range_y') 
        xlabel = Dropdown(
                  description='* X-Axis',
                  values=controller.ds.independent_variables, 
                  options=controller.ds.independent_variables, 
                  value=xaxis)
        ylabel = Dropdown(
                  description='* Y-Axis',
                  values=['None'] + controller.ds.independent_variables,
                  options=['None'] + controller.ds.independent_variables,
                  value=yaxis)
        xmin = FloatText(description='* X-Min',
                                       value=range_x[0])
        xmax = FloatText(description='* X-Max',
                                       value=range_x[1])
        ymin = FloatText(description='* Y-Min',
                                       value=range_y[0])
        ymax = FloatText(description='* Y-Max',
                                       value=range_y[1])
        center_axes = Checkbox(description='Center Axes',
                                             value=center)
        boundaries = Checkbox(description='Draw Boundaries',
                                            value=False)
        plot_type = Dropdown(description='* Plot Type',
                             values=self.widget_types,
                             options=self.widget_types,
                             value='Design Space (interactive)')
        title_widget = Text(description='Title')
        caption_widget = Textarea(description='Caption')
        included = controller.defaults('included_cases')
        if included is None:
            included = []
        if isinstance(included, list) is False:
            included = [included]
        included_widget = Textarea(description='Only Cases',
                                                 value=', '.join(included))
        wi = VBox(children=[xlabel, ylabel, plot_type, 
                                               xmin, xmax, ymin, ymax,
                                               center_axes, boundaries,
                                               title_widget, caption_widget,
                                               included_widget])
        for i in [xlabel, ylabel, plot_type]:
            i.on_trait_change(self.update_field, 'value')    
        plot_type.widget_container = wi
        button = Button(value=False, description='Add Plot')

        button.on_click(self.make_plot)
        button.xlabel = xlabel
        button.ylabel = ylabel
        button.xmin = xmin
        button.xmax = xmax
        button.ymin = ymin
        button.ymax = ymax
        button.center_axes = center_axes
        button.boundaries = boundaries
        button.plot_type = plot_type
        button.title = title_widget
        button.caption = caption_widget
        button.included = included_widget
        button.wi = wi
        self.show_colorbar = Checkbox(description='Split Colorbar', value=False)
        button.show_colorbar = self.show_colorbar
        self.title = title_widget
        self.caption = caption_widget
        self.boundaries = boundaries
        self.plot_type = plot_type
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.ymin = ymin
        self.ymax = ymax
        self.xmin = xmin
        self.xmax = xmax
        self.button_make_interactive = button # Added by Miguel
        
        self.button_make_interactive.upper_title = 'Design Space Plot'
        self.button_make_interactive.colorbar = False
        self.button_make_interactive.highlight_flag = False

        box_add_plot = Box()
        box_add_plot.children = [button]
        
        add_plot = VBox(description='Add Plot',
                                           children=[wi,
                                                     self.plot_data, 
                                                     box_add_plot])
        self.update_plot_widget('value', 'Design Space (Interactive)')
        return ('Create Plot', add_plot)
        
    def update_field(self, name, value):
        controller = self.controller
        self.is_1D = True if str(self.ylabel.value) == 'None' else False
        controller.set_defaults('xaxis', str(self.xlabel.value))
        controller.set_defaults('yaxis', str(self.ylabel.value))
        controller.set_defaults('range_x', [self.xmin.value, self.xmax.value])
        controller.set_defaults('range_y', [self.ymin.value, self.ymax.value])
        if self.is_1D is True:
            self.ymin.visible = False
            self.ymax.visible = False
            self.boundaries.visible = False
            self.ylabel.description = '* 2D plot'
        else:
            self.ymin.visible = True
            self.ymax.visible = True
            self.boundaries.visible = True
            self.ylabel.description = '* Y-Axis'
        self.update_plot_widget(name, value)
        
    def stability_2D_plot_widget(self):
        controller = self.controller
        resolution_widget = FloatText(description='Resolution', value=controller.defaults('resolution'))
        colorbar_widget = Checkbox(description='Show colorbar', value=True)
        wi = VBox(children=[resolution_widget, colorbar_widget])
        wi.resolution = resolution_widget
        wi.colorbar_widget = colorbar_widget
        self.plot_data.children = [wi]
        self.title.value = 'System design space showing stability of the fixed points'
        self.caption.value = 'Number of eigenvalues with positive real part represented as a heat map on the z-axis.'
        return
    
    def stability_1D_plot_widget(self):
        controller = self.controller
        stability_1D_options = ['log('+ i + ')' for i in controller.ds.dependent_variables]
        zlim = controller.defaults('zlim')
        function_widget = Dropdown(description='* Y-Axis',
                               values=stability_1D_options, value=stability_1D_options[0]) #Modified by Miguel to fix drop-down menu. #value = 'log('+controller.ds.dependent_variables[0]+')')
        resolution_widget = FloatText(description='Resolution', value=controller.defaults('resolution'))
        zlim_auto = (zlim is None)
        zlim_widget = Checkbox(description='Automatic Y-Lim', value=zlim_auto)
        if zlim_auto is True:
            zlim = [0., 0.]
        zmin_widget = FloatText(description='Y-Min', value=zlim[0])
        zmax_widget = FloatText(description='Y-Max', value=zlim[1])
        wi = VBox(children=[function_widget, resolution_widget,
                                               zlim_widget, zmin_widget, zmax_widget,
                                               ])
        wi.function = function_widget
        wi.resolution = resolution_widget
        wi.zlim = zlim_widget
        wi.zmin = zmin_widget
        wi.zmax = zmax_widget
        self.plot_data.children = [wi]
        self.title.value = 'System design space showing stability of the fixed points'
        self.caption.value = 'Number of eigenvalues with positive real part represented by line style: '
        self.caption.value += '0 eigenvalues w/ positive real part (solid); '
        self.caption.value += '1 eigenvalues w/ positive real part  (red dotted); '
        self.caption.value += '2 eigenvalues w/ positive real part (green dashed).'
        return
        
    def eigenvalue_2D_plot_widget(self):
        controller = self.controller
        zlim = controller.defaults('zlim')
        component_widget = Dropdown(description='Complex component',
                                    values=['Real', 'Imaginary'],
                                    options=['Real', 'Imaginary'],
                                    value='Real')
        resolution_widget = FloatText(description='Resolution', value=controller.defaults('resolution'))
        zlim_auto = (zlim is None)
        zlim_widget = Checkbox(description='Automatic Z-Lim', value=zlim_auto)
        if zlim_auto is True:
            zlim = [0., 0.]
        zmin_widget = FloatText(description='Z-Min', value=zlim[0])
        zmax_widget = FloatText(description='Z-Max', value=zlim[1])
        parallel_widget = Checkbox(description='Compute in Parallel', value=False)
        number_dynamic = len(controller.ds.dependent_variables)
        number_dynamic -= len(controller.ds.auxiliary_variables)
        select_widget = Dropdown(
                         description='Rank to Plot',
                         values = [str(i+1) for i in range(number_dynamic)],
                         options = [str(i+1) for i in range(number_dynamic)],
                         value=str(1))
        wi = VBox(children=[component_widget, 
                                               select_widget,
                                               resolution_widget,
                                               zlim_widget,
                                               zmin_widget,
                                               zmax_widget,
                                               parallel_widget])
        wi.component = component_widget
        wi.select = select_widget
        wi.resolution = resolution_widget
        wi.parallel = parallel_widget
        wi.zlim = zlim_widget
        wi.zmin = zmin_widget
        wi.zmax = zmax_widget
        self.plot_data.children = [wi]
        self.title.value = 'System design space showing the dominant eigenvalue of the fixed points'
        self.caption.value = 'Dominant eigenvalue represented as a heat map on the z-axis.'
        return
    
    def function_2D_plot_widget(self):
        controller = self.controller
        zlim = controller.defaults('zlim')
        value = str(self.plot_type.value)
        log_linear_widget = Checkbox(description='Function is log linear',
                                                   value=True)
        if value == 'Steady State Flux':
            flux_options = ['log(V_'+ i + ')' for i in controller.ds.dependent_variables_no_algebraic_constraints]
            function_widget = Dropdown(values=flux_options,
                                       options=flux_options)
            self.title.value = 'System design space showing a steady state flux'
            self.caption.value = 'Steady state flux shown as a heat map on the z-axis.'
        elif value == 'Steady State Function':
            function_widget = Text(description='Function', value='')
            log_linear_widget.value = False
            self.title.value = 'System design space showing a function at steady state'
            self.caption.value = 'Steady state function shown as a heat map on the z-axis.'
        else:
            ss_options = ['log('+ i + ')' for i in controller.ds.dependent_variables_no_algebraic_constraints]
            function_widget = Dropdown(values=ss_options,
                                       options=ss_options)
            self.title.value = 'System Design Space showing a steady state concentration'
            self.caption.value = 'Steady state concentration shown as a heat map on the z-axis.'
        resolution_widget = FloatText(description='Resolution', value=controller.defaults('resolution'))
        parallel_widget = Checkbox(description='Compute in Parallel', value=False)
        zlim_auto = (zlim is None)
        zlim_widget = Checkbox(description='Automatic Z-Lim', value=zlim_auto)
        if zlim_auto is True:
            zlim = [0., 0.]
        zmin_widget = FloatText(description='Z-Min', value=zlim[0])
        zmax_widget = FloatText(description='Z-Max', value=zlim[1])
        wi = VBox(children=[function_widget, resolution_widget,
                                               zlim_widget, zmin_widget, zmax_widget,
                                               parallel_widget, log_linear_widget])
        wi.function = function_widget
        wi.resolution = resolution_widget
        wi.parallel = parallel_widget
        wi.log_linear = log_linear_widget
        wi.zlim = zlim_widget
        wi.zmin = zmin_widget
        wi.zmax = zmax_widget
        self.plot_data.children = [wi]
        return

        
    def function_1D_plot_widget(self):
        controller = self.controller
        zlim = controller.defaults('zlim')
        value = str(self.plot_type.value)
        self.function_2D_plot_widget()
        wi = self.plot_data.children[0]
        if value == 'Steady State Flux':
            self.title.value = 'System design space showing a steady state flux'
            self.caption.value = 'Steady state flux shown on the y-axis.'
        elif value == 'Steady State Function':
            self.title.value = 'System design space showing a function at steady state'
            self.caption.value = 'Steady state function shown on the y-axis.'
        else:
            self.title.value = 'System Design Space showing a steady state concentration'
            self.caption.value = 'Steady state concentration shown on the y-axis.'
        wi.function.description = 'Y-Axis'
        wi.zlim.description = 'Automatic Y-Lim'
        wi.zmin.description = 'Y-Min'
        wi.zmax.description = 'Y-Max'
        return
            
    def update_plot_widget(self, name, value):
        controller = self.controller
        zlim = controller.defaults('zlim')
        value = str(self.plot_type.value)
        if value == 'Design Space (interactive)':
            wi = VBox(children=[self.show_colorbar])
            self.plot_data.children = [wi]
            self.title.value = 'System design space'
            self.caption.value = 'System design space with the enumerated qualitatively-distinct phenotypes represented on the z-axis, identified by color.'                
        elif value == 'Design Space':
            intersections_widget = Dropdown(description='# Intersections',
                                            values=['Single',
                                                    'Single and Triple',
                                                    'Triple',
                                                    'All',],
                                            options=['Single',
                                                     'Single and Triple',
                                                     'Triple',
                                                     'All',],
                                            value='Single')
            
            ### Add Advanced Settings buttons. Create a container and place all buttons here.
            advanced_button = Button(description='Open Advanced Settings')
            advanced_button.on_click(self.expand_advanced_settings)
            highlight_phenotypes = Text(description = 'Highlight Phenotype:', visible = False)
            figure_title = Text(description= 'Upper Figure Title', value='Design Space Plot', visible = False)
            show_colorbar = Checkbox(description='Show color bar', value=False, visible = False)
            accept_advanced = Button(description='Accept Changes', visible = False)
            advanced_container = Box(children = [highlight_phenotypes, figure_title, show_colorbar, accept_advanced, advanced_button])
            # Attach widgets to advanced_button
            advanced_button.phenotypes = highlight_phenotypes
            advanced_button.title = figure_title
            advanced_button.show_colorbar = show_colorbar
            advanced_button.accept = accept_advanced
            advanced_button.state = 'open'
            ## Attach widgets to accept_advanced button and define on_click action.
            accept_advanced.phenotypes = highlight_phenotypes
            accept_advanced.title = figure_title
            accept_advanced.colorbar = show_colorbar
            accept_advanced.advanced_button = advanced_button
            accept_advanced.on_click(self.update_static_variables)
            
            wi = VBox(children=[intersections_widget, advanced_container])
            wi.intersections = intersections_widget
            self.title.value = 'System design space'
            self.caption.value = 'Enumerated qualitatively-distinct phenotypes represented on the z-axis and identified by color.'                
            self.plot_data.children = [wi]
        elif value == 'Stability':
            if self.is_1D:
                self.stability_1D_plot_widget()
            else:
                self.stability_2D_plot_widget()

        elif value == 'Eigenvalues':
            if self.is_1D:
                return
            else:
                self.eigenvalue_2D_plot_widget()
        else:
            if self.is_1D:
                self.function_1D_plot_widget()
            else:
                self.function_2D_plot_widget()
        if controller.name != '':
            title = 'Analysis of the ' + controller.name + ' by ' + self.title.value.lower()
            self.title.value = title
        
        
    def expand_advanced_settings(self, b): # Miguel

        if b.state == 'open':
            b.description = 'Cancel'
            b.phenotypes.visible = True
            b.title.visible = True
            b.show_colorbar.visible = True
            b.accept.visible = True
            b.state = 'close'
        else:
            b.description = 'Advanced Settings'
            b.visible = True
            b.phenotypes.visible = False
            b.title.visible = False
            b.show_colorbar.visible = False
            b.accept.visible = False
            b.state = 'open'

    def update_static_variables(self,b):
        # close advanced settings menu
        b.phenotypes.visible = False
        b.included = b.phenotypes
        b.title.visible = False
        b.colorbar.visible = False
        b.visible = False
        b.advanced_button.description = 'Open Advanced Settings'
        b.advanced_button.state = 'open'
        
        # Update Fields
        self.button_make_interactive.upper_title = b.title.value
        # show color bar
        self.button_make_interactive.colorbar = b.colorbar.value
        # Hihlight cases
        self.button_make_interactive.highlight_cases = b.phenotypes.value
        if str(self.button_make_interactive.highlight_cases) != '':
            self.button_make_interactive.highlight_flag = True
        else:
            self.button_make_interactive.highlight_flag = False

    def make_plot(self, b):
        controller = self.controller
        b.description = 'Creating plot... Please Wait.'
        b.disabled = True
        # if 1 == 1:
        try:
            b.pvals = controller.pvals.copy()
            if b.plot_type.value == 'Design Space (interactive)':
                self.make_interactive_plot(b)
            elif b.plot_type.value == 'Design Space':
                self.make_static_plot(b)
            elif b.plot_type.value == 'Stability':
                self.make_stability_plot(b)
            elif b.plot_type.value == 'Eigenvalues':
                self.make_eigenvalue_plot(b)
            else:
                self.make_function_plot(b)
        except Exception as e:
            print(e)    # print exception for debugging purposes.
            close_button = Button(description="Close")
            error_message = '<div width="200px" style="float:center; text-align:center;">'
            error_message += '<b>An error occured while plotting'
            if str(e).count('Case') == 1:
                error_message += ': ' + str(e) + '</b></div>'
            else:
                error_message += '</b></div>'
            error_window = Popup(children=[HTML(value=error_message),close_button])
            close_button.window = error_window
            close_button.on_click(lambda x: x.window.close())
            if old_ipython is False:
                error_window.box_style='danger'
                close_button.float = 'center'
                error_window.width='250px'
                error_window.height='150px'
            display(error_window)
            b.description = 'Add Plot'
            b.disabled = False
        finally:
            b.description = 'Add Plot'
            b.disabled = False

    def axes_ranges(self, b):
        pvals = self.controller.pvals
        ranges = [[b.xmin.value, b.xmax.value],[b.ymin.value, b.ymax.value]]
        if b.center_axes.value is False:
            return ranges
        if str(b.ylabel.value) != 'None':
            ranges = [[pvals[str(b.xlabel.value)]*i for i in ranges[0]],
                      [pvals[str(b.ylabel.value)]*i for i in ranges[1]]]
        else:
            ranges = [[pvals[str(b.xlabel.value)]*i for i in ranges[0]],
                      None]

        return ranges
    
    def included_cases(self, b):
        included = str(b.included.value).strip()
        if len(included) == 0:
            return None
        included = [i.strip() for i in included.split(',')]
        return included
        
    def make_interactive_plot(self, b):
        controller = self.controller
        if str(b.ylabel.value) == 'None':
            return
        button = Button(description='Stop interactive plot')
        button2 = Button(description='Adjust Sliders')
        self.sliders_menu = Box()
        self.sliders_menu.children = [button2]
        button2.on_click(self.customize_sliders)
        button.on_click(self.remove_plot)
        button.name = 'Interactive Plot (' + str(np.random.randint(0, 1000)) + ')'
        self.button_make_interactive.name = button.name

        if self.show_colorbar.value is True:
            image_widget_fig = Image()
            image_widget_color = HTML()
            image_widget_fig.description = 'Figure'
            image_widget_fig.width = '100%'
            image_widget_color.description = 'Colorbar'
            image_widget = Tab(children=[VBox(children=[image_widget_fig]), image_widget_color])
            image_widget.set_title(0, 'Figure')
            image_widget.set_title(1, 'Colorbar')
            image_widget.selected_index = 0
            image_widget.description = 'Figure_No_Colorbar'
        else:
            image_widget = Image()
            image_widget.width = '100%'
            image_widget.description = 'Figure'

        popup_widget = Popup(children=[image_widget])
        rangex, rangey = self.axes_ranges(b)
        
        # initialize dictionary
        slider_dictionary = {i: [1e-5, 1e5] for i in controller.ds.independent_variables}
        
        # Update values for x
        #slider_dictionary.update({b.xlabel.value:[b.xmin.value, b.xmax.value]})

        # update values for y
        #slider_dictionary.update({b.ylabel.value:[b.ymin.value, b.ymax.value]})
        
        # if sliders were modified, update dictionary.
        if hasattr(b, "modify_sliders") and b.modify_sliders is True:
            #update slider dictionary.
            #slider_dictionary.update({b.sliders_variable.value:[10**b.sliders_min.value , 10**b.sliders_max.value]})
            slider_dictionary.update(b.sliders_update_dictionary)
            
        interactive_plot = controller.ds.draw_2D_slice_notebook(controller.pvals, str(b.xlabel.value),
                                                                str(b.ylabel.value),
                                                                rangex, rangey,
                                                                slider_dictionary,
                                                                intersections=[1, 3], #[1, 3]
                                                                image_container=image_widget)
        wi = VBox(description=button.name,
                  children=[interactive_plot,
                            button, self.sliders_menu,
                            popup_widget])
        controller.options.update({'xaxis': str(b.xlabel.value),
                                   'yaxis': str(b.ylabel.value),
                                   'x_range': rangex,
                                   'y_range': rangey})
        controller.update_child(button.name, wi)
        
    def make_static_plot(self, b):
        if hasattr(b, 'update'):
            b = b.b
        controller = self.controller
        if old_ipython is True:
            if self.button_make_interactive.colorbar is False:
                fig = plt.figure(figsize=[5.5698, 4], dpi=600, facecolor='w') # original size of figure is 7,4
            else:
                fig = plt.figure(figsize=[7, 4], dpi=600, facecolor='w')
            
            ax = fig.add_axes([0.1714, 0.2, 0.6, 0.7])
            ax.set_title(b.upper_title) # 'Design Space plot'
        else:
            fig = plt.figure(figsize=[5, 4], dpi=600, facecolor='w')
            ax = fig.add_axes([0.24, 0.2, 0.72, 0.7])
            ax.set_title('Design Space plot')
        plot_data = self.plot_data.children[0]
        intersects = plot_data.intersections.value
        intersections_dict = {'Single':[1],
                              'Single and Triple': [1, 3],
                              'Triple':[3],
                              'All':range(1, 100)}
        rangex, rangey = self.axes_ranges(b)
        ec = 'k' if b.boundaries.value is True else 'none'
        if str(b.ylabel.value) != 'None':
            colorbar = self.button_make_interactive.colorbar if old_ipython is True else False
            colors=controller.ds.draw_2D_slice(ax, controller.pvals,
                                               str(b.xlabel.value), str(b.ylabel.value),
                                               rangex, rangey,
                                               intersections=intersections_dict[intersects],
                                               included_cases=self.included_cases(b),
                                               ec=ec,
                                               colorbar=colorbar) # changed from colorbar=colorbar to colorbar = False
            if self.button_make_interactive.highlight_flag:
                case = controller.ds(str(self.button_make_interactive.highlight_cases))
                case.draw_2D_slice(ax, controller.pvals, str(b.xlabel.value), str(b.ylabel.value),
                                       rangex, rangey, fc='none', ec='k', lw='2.')
        else:
            colors=controller.ds.draw_1D_slice(ax, controller.pvals, str(b.xlabel.value),
                                               rangex,
                                               intersections=intersections_dict[intersects],
                                               included_cases=self.included_cases(b))
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, 
                                      title=b.title.value,
                                      caption=b.caption.value,
                                      pvals=b.pvals,
                                      colors=colors)
        controller.options.update({'xaxis':str(b.xlabel.value),
                                   'yaxis':str(b.ylabel.value),
                                   'x_range':rangex, 
                                   'y_range':rangey,
                                   'included_cases':self.included_cases(b)})
        plt.close()
        
    def make_stability_plot(self, b):
        controller = self.controller
        if hasattr(self.plot_data.children[0], 'colorbar_widget'):
            show_colorbar = self.plot_data.children[0].colorbar_widget.value
        else:
            show_colorbar = True
        if show_colorbar is True:
            fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
        else:
            fig = plt.figure(figsize=[6*0.8, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        ax.set_title('Stability plot')
        plot_data = self.plot_data.children[0]
        resolution = plot_data.resolution.value
        rangex, rangey = self.axes_ranges(b)
        ec = 'k' if b.boundaries.value is True else 'none'
        if str(b.ylabel.value) != 'None':
            cf, colors = controller.ds.draw_2D_positive_roots(ax, controller.pvals,
                                                              str(b.xlabel.value), str(b.ylabel.value),
                                                              rangex, rangey,
                                                              colorbar=show_colorbar,
                                                              resolution=resolution,
                                                              included_cases=self.included_cases(b)
                                                             )
            if ec == 'k':
                controller.ds.draw_2D_slice(ax, controller.pvals, str(b.xlabel.value), str(b.ylabel.value),
                                            rangex, rangey,
                                            intersections=[1],
                                            included_cases=self.included_cases(b),
                                            colorbar=False,
                                            facecolor='none',
                                            ec=ec)
        else:
            zlim = None
            function = str(plot_data.function.value)
            if plot_data.zlim.value == False:
                zlim = [plot_data.zmin.value, plot_data.zmax.value]
            controller.ds.draw_1D_positive_roots(ax, function, controller.pvals, 
                                                 str(b.xlabel.value), rangex,
                                                 ylim=zlim,
                                                 resolution=resolution)
            controller.set_defaults('zlim', zlim)
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()

        if show_colorbar is False and str(b.ylabel.value) != 'None':
            controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value, pvals=b.pvals,
                                          colors=colors)
        else:
            controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value, pvals=b.pvals)

        plt.close()
        controller.options.update({'xaxis':str(b.xlabel.value),
                                   'yaxis':str(b.ylabel.value),
                                   'x_range':rangex, 
                                   'y_range':rangey,
                                   'included_cases':self.included_cases(b),
                                   'resolution':resolution})
        
    
    def make_function_plot(self, b):
        controller = self.controller
        plot_data = self.plot_data.children[0]
        log_linear = plot_data.log_linear.value
        function = str(plot_data.function.value)
        resolution = plot_data.resolution.value
        parallel = plot_data.parallel.value
        zlim = None
        if plot_data.zlim.value == False:
            zlim = [plot_data.zmin.value, plot_data.zmax.value]
        fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        fn = dspace.Expression(function)
        rangex, rangey = self.axes_ranges(b)
        ax.set_title('$' + fn.latex(substitution_dictionary=controller.symbols) + '$')
        ec = 'k' if b.boundaries.value is True else 'none'
        if str(b.ylabel.value) != 'None':
            controller.ds.draw_2D_ss_function(ax, function, controller.pvals, 
                                              str(b.xlabel.value),
                                              str(b.ylabel.value),
                                              rangex, rangey, zlim=zlim,
                                              log_linear=log_linear, resolution=resolution, 
                                              parallel=parallel,
                                              included_cases=self.included_cases(b))
            if ec == 'k':
                controller.ds.draw_2D_slice(ax, controller.pvals, str(b.xlabel.value), str(b.ylabel.value),
                                            rangex, rangey,
                                            intersections=[1],
                                            included_cases=self.included_cases(b),
                                            colorbar=False,
                                            facecolor='none',
                                            ec=ec)
        else:
            controller.ds.draw_1D_ss_function(ax, function, controller.pvals, 
                                              str(b.xlabel.value),
                                              rangex, ylim=zlim,
                                              resolution=resolution, 
                                              included_cases=self.included_cases(b))
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value, pvals=b.pvals)
        plt.close()
        controller.options.update({'xaxis':str(b.xlabel.value),
                                   'yaxis':str(b.ylabel.value),
                                   'x_range':rangex, 
                                   'y_range':rangey,
                                   'included_cases':self.included_cases(b),
                                   'resolution':resolution,
                                   'zlim':zlim})
        
    def make_eigenvalue_plot(self, b):
        controller = self.controller
        if str(b.ylabel.value) == 'None':
            return
        plot_data = self.plot_data.children[0]
        component = str(plot_data.component.value)
        resolution = plot_data.resolution.value
        parallel = plot_data.parallel.value
        rank = int(plot_data.select.value)
        zlim = None
        if plot_data.zlim.value == False:
            zlim = [plot_data.zmin.value, plot_data.zmax.value]
        fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        rangex, rangey = self.axes_ranges(b)
        ax.set_title('Dominant Eigenvalue ('+component+')')
        controller.ds.draw_2D_dominant_eigenvalues(ax, controller.pvals, 
                                                   str(b.xlabel.value),
                                                   str(b.ylabel.value),
                                                   rangex, rangey, zlim=zlim,
                                                   component=component.lower(),
                                                   resolution=resolution, 
                                                   parallel=parallel,
                                                   included_cases=self.included_cases(b),
                                                   cmp=lambda eig : 
                                                       eigenvalue_compare(eig,
                                                                          component=component.lower(),
                                                                          rank=rank))
        ec = 'k' if b.boundaries.value is True else 'none'
        if ec == 'k':
            controller.ds.draw_2D_slice(ax, controller.pvals, str(b.xlabel.value), str(b.ylabel.value),
                                        rangex, rangey,
                                        intersections=[1],
                                        included_cases=self.included_cases(b),
                                        colorbar=False,
                                        facecolor='none',
                                        ec=ec)
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value, pvals=b.pvals)
        plt.close()
        controller.options.update({'xaxis':str(b.xlabel.value),
                                   'yaxis':str(b.ylabel.value),
                                   'x_range':rangex, 
                                   'y_range':rangey,
                                   'included_cases':self.included_cases(b),
                                   'resolution':resolution,
                                   'zlim':zlim})
                
    def remove_plot(self, b):
        controller = self.controller
        controller.update_child(b.name, None)

    def customize_sliders(self, b):
        controller = self.controller
        selectVariable = Dropdown(description = 'Select Variable', values=controller.ds.independent_variables, options=controller.ds.independent_variables)  ## Figure out how to sort this!
        selectMin = FloatText(description='Min value (log scale)', value=-5)
        selectMax = FloatText(description='Max value (log scale)', value=5)
        selectStep = FloatText(description='Step (log scale)', value=0.5)

        accept_button = Button(description='Accept')
        accept_button.selectMin = selectMin
        accept_button.selectMax = selectMax
        accept_button.selectStep = selectStep
        accept_button.variable = selectVariable
        accept_button.on_click(self.adjustslider)
        
        cancel_button = Button(description='Cancel')
        cancel_button.on_click(self.redoslidersbutton)
        
        self.sliders_menu.children = [selectVariable,selectMin, selectMax, accept_button, cancel_button] # step option was removed. add selectStep to consider it again.

    def redoslidersbutton(self,b):
    
        button2 = Button(description ='Adjust Sliders')
        self.sliders_menu.children = [button2]
        button2.on_click(self.customize_sliders)

    def adjustslider(self,b):  #Tag
        #Stop interactive plot
        self.remove_plot(self.button_make_interactive)
        
        # Extract values from button.
        min_value = b.selectMin
        max_value = b.selectMax
        variable = b.variable
        
        # calculate step if not provided
        if b.selectStep.value == 0:
            b.selectStep.value = (float(max_value) - float(min_value))/20
            step = b.selectStep.value
        else:
            step = b.selectStep.value
        
        # pass slider parameters
        self.button_make_interactive.sliders_min = min_value
        self.button_make_interactive.sliders_max = max_value
        self.button_make_interactive.sliders_step = step
        self.button_make_interactive.modify_sliders = True
        self.button_make_interactive.sliders_variable = b.variable
        
        if not hasattr(self.button_make_interactive, "sliders_update_dictionary"):
            setattr(self.button_make_interactive, 'sliders_update_dictionary',{})
        self.button_make_interactive.sliders_update_dictionary.update({variable.value: [10**min_value.value, 10**max_value.value]})
        
        # Generate Widget
        self.make_interactive_plot(self.button_make_interactive)



class DisplayFigures(object):
                         
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'figures_widget', None)
        setattr(self, 'unsaved', None)
        setattr(self, 'delete_menu', None)
        
    def create_figures_widget(self):
        
        controller = self.controller
        self.figures_widget = VBox()
        self.unsaved = VBox()
        self.delete_figures_button = Button (description = "Delete Selected Figures")   #Miguel
        self.delete_menu = VBox(children = [self.delete_figures_button])                #Miguel
        self.delete_figures_button.on_click(self.open_delete_menu)                      #Miguel
        unsaved = '<center><b style="color:red;">Figures that will not be saved:</b></center><br><hr>'
        self.figures = VBox(children=[self.delete_menu, self.figures_widget, #Miguel
                                                         HTML(value=unsaved),
                                                         self.unsaved])
        controller.update_child('Figures', self.figures)
        
    def add_figure(self, image_data, title='', caption = '', pvals=None, colors=None):
        
        controller = self.controller
        
        if pvals is not None and len(pvals) <= 40:
            caption += ' Figure generated with the following parameter values: '
#            # if the length of the dictionary pvals is different from controller.pvals, list the elements of pvals, otherwise, list elements of controller.pvals.
#            if len(pvals) != len(controller.pvals):
#                caption += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals.keys())])
#            else:
            caption += '; '.join([i + ' = ' + str(controller.pvals[i]) for i in sorted(controller.pvals.keys())])
            # expand caption here. Add kinetic order and constraints.
            if str(controller.replacements_caption.value) != '':
                caption += ';' + ' Kinetic order(s): ' + str(controller.replacements_caption.value)
            if str(controller.constraints_caption.value) != '':
                caption += ';' + ' Parametric constraints: ' + str(controller.constraints_caption.value)
            caption += '.'
            if hasattr(controller, 'initial_conditions_caption'):
                if controller.initial_conditions_caption != '':
                    caption += ' Initial conditions are: ' + controller.initial_conditions_caption + '.'
                    controller.initial_conditions_caption = ''

        self.add_figure_widget(image_data, title=title, caption=caption, pvals=pvals, colors=colors)
        
    def remove_unsaved_figure(self, b):
        children = [i for i in self.unsaved.children] 
        children.remove(b.wi)
        self.unsaved.children = children
        
    def save_unsaved_figure(self, b):
        controller = self.controller
        self.remove_unsaved_figure(b)        
        self.save_figure(b.image_data, title=b.title, caption=b.caption, pvals = b.pvals, colors=b.colors)
        controller.save_widget_data(b)
        
    def save_figure(self, image_data, title='', caption = '', pvals=None, colors=None):
        controller = self.controller
        figures = controller.figure_data
        figures.append((image_data, title, caption, pvals, colors))
        self.save_figure_widget(image_data, title=title, 
                                caption=caption, pvals=pvals, colors=colors)
        
    def colorbar_tabs_html(self, colors, height=None):
        tab_dicts = {}
        html_widgets = []
        for i in colors:
            key=len(i.split(','))
            if key not in tab_dicts:
                tab_dicts[key] = {}
            tab_dicts[key][i]  = '#%02x%02x%02x' % tuple([j*255 for j in colors[i][:3]])
        keys = sorted(tab_dicts)
        labels = [sorted(tab_dicts[i]) for i in keys]
        lengths = [len(tab_dicts[i]) for i in keys]
        max_length = max(lengths)
        html_str = '<table style="border:0;"' + ('height="'+height+'">' if height is not None else '>')
        for i in range(max_length):
            html_str += '<tr style="border:0;">'
            for j in range(len(labels)):
                if i < lengths[j]:
                    key = keys[j]
                    label = labels[j][i]
                    html_str += '<td style="border:0; background-color:{0}; padding:10px;" />'.format(tab_dicts[key][label])
                    html_str += '<td style="border:0; white-space: nowrap; font-size:100%">'+label+'</td>'
                else:
                    html_str += '<td style="border:0;" />'
                    html_str += '<td style="border:0;" />'
            html_str += '</tr>'
        html_str += '</table>'
        return html_str
    
    def colorbar_tabs(self, colors):
        html_str = self.colorbar_tabs_html(colors)
        html_ob = HTML(value=html_str)
        html_ob.set_css('height', '400px')
        html_widgets = [html_ob]
        return html_widgets

    def _image_widget_new(self, image_data, colors=None):
        b64_data = base64.b64encode(image_data)
        if colors is None:
            html_str = '<img src="data:image/png;base64,'+b64_data+'" width="100%" />'
        else:
            html_str = '''
<style>.col1 {display: none; }
.col2 {display: none; }
.col3 {display: none; }
table.show1 .col1 { display: table-cell; }
table.show2 .col2 { display: table-cell; }
table.show3 .col3 { display: table-cell; }
</style>
<script>
function toggleColumn(el, n) {
    var currentClass = el.parentElement.parentElement.parentElement.parentElement.className;
    if (currentClass.indexOf("show"+n) != -1) {
        el.parentElement.parentElement.parentElement.parentElement.className = currentClass.replace("show"+n, "");
    }
    else {
        el.parentElement.parentElement.parentElement.parentElement.className += " " + "show"+n;
    }
}
</script>
'''
            html_str += '<table style="border:0;" width="100%"><tr>'
            html_str += '<th><button onclick="toggleColumn(this,2)">Toggle Colorbar</button></th><th></th></tr>'
            html_str += '<tr style="border:0;paddding-right:5px;"><td style="border:0;">'
            html_str += '<img src="data:image/png;base64,'+b64_data+'" / width="100%"></td>'
            html_str += '<td style="border:0;" class="col2">'+self.colorbar_tabs_html(colors)+'</td></tr></table>'
        image_widget = HTML(value=html_str)
        image_widget.width='100%'
        return image_widget
        
    def add_figure_widget(self, image_data, title='', caption = '', pvals=None, colors=None):
        children = [i for i in self.unsaved.children]      
        if len(title) > 0:
            title = title + '.'
        if len(caption) > 0:
            caption = '  ' + caption
        html_str = '<b>'+title+'</b>' + caption
        html_widget = HTML(value=html_str)
        save_button = Button(description='Save Figure')
        save_button.image_data = image_data
        save_button.title = title
        save_button.caption = caption
        save_button.on_click(self.save_unsaved_figure)
        save_button.pvals = pvals
        save_button.colors = colors
        close_button = Button(description='Remove Figure')
        close_button.on_click(self.remove_unsaved_figure)
        restore_pvals = Button(description='Restore Parameter Values')
        restore_pvals.pvals = pvals
        if pvals is None:
            restore_pvals.visible = False
        if old_ipython is True:
            image_widget = Image()
            image_widget.value = image_data
            tab_widget = VBox(children=[image_widget, html_widget])
            if colors is not None:
                html_widgets = self.colorbar_tabs(colors)
                tab_widget.description='Figure'
                tab_widget = Tab(children=[tab_widget]+html_widgets)
                # tab_widget.set_css('height', '500px')
        else:
            image_widget = self._image_widget_new(image_data, colors=colors)
            tab_widget = VBox(children=[image_widget, html_widget])
            toggle = Button(description='Hide')
            
        restore_pvals.on_click(self.restore_figure_pvals)
        
        if old_ipython is True:
            wi = Popup(children=[close_button, save_button, tab_widget, restore_pvals])
        else:
            contents = VBox(children=[close_button, save_button, tab_widget, restore_pvals])
            contents.width = '100%'
            wi = Popup(children=[contents])
            wi.border = '1px solid black'
        save_button.wi = wi
        close_button.wi = wi
        children.append(wi)
        self.unsaved.children = children
        if colors is not None:
            if old_ipython is True:
                tab_widget.set_title(0, 'Figure')
                tab_widget.set_title(1, 'Colorbar')
            
    def save_figure_widget(self, image_data, title='', caption = '', pvals=None, colors=None):
        image_widget = Image()
        image_widget.value = image_data
        children = [i for i in self.figures_widget.children]      
        html_str = '<b>Figure '+str(len(children)+1)+'.  '+title+'</b>' + caption
        html_widget = HTML(value=html_str)
        restore_pvals = Button(description='Restore Parameter Values')
        restore_pvals.pvals = pvals
        if pvals is None:
            restore_pvals.visible = False
        if old_ipython is True:
            image_widget = Image()
            image_widget.value = image_data
            tab_widget = VBox(children=[image_widget, html_widget])
            if colors is not None:
                html_widgets = self.colorbar_tabs(colors)
                tab_widget.description='Figure'
                tab_widget = Tab(children=[tab_widget]+html_widgets)
        else:
            image_widget = self._image_widget_new(image_data, colors=colors)
            tab_widget = VBox(children=[image_widget, html_widget])
        restore_pvals.on_click(self.restore_figure_pvals)
        if old_ipython is True:
            wi = Popup(children=[tab_widget, restore_pvals])
        else:
            contents = VBox(children=[tab_widget, restore_pvals])
            contents.width = '100%'
            wi = Popup(children=[contents])
            wi.border = '1px solid black'
        children.append(wi)
        self.figures_widget.children = children
        if colors is not None:
            if old_ipython is True:
                tab_widget.set_title(0, 'Figure')
                tab_widget.set_title(1, 'Colorbar')
                
    def restore_figure_pvals(self, b):
        controller = self.controller
        controller.pvals = b.pvals

    def load_widgets(self):
        controller = self.controller
        for data in controller.figure_data:
            if len(data) == 3:
                self.save_figure_widget(data[0], title=data[1], caption=data[2])
            elif len(data) == 4:
                self.save_figure_widget(data[0], title=data[1], caption=data[2], pvals=data[3])
            elif len(data) == 5:
                self.save_figure_widget(data[0], title=data[1], caption=data[2], pvals=data[3], colors=data[4])

    def open_delete_menu(self, b): #Miguel
        description = Latex()
        description.value = "Specify figures to delete. Please note that deleted items cannot be recovered!"
        indText = Textarea(description="Individual figures (comma separated values)")
        rangeText = Textarea(description="Range of figures (e.g 1-2)")
        deleteButton = Button(description="Delete figures")
        cancelButton = Button(description="Cancel")
        self.delete_menu.children = [description, indText, rangeText, deleteButton, cancelButton]
        deleteButton.indValue = indText
        deleteButton.rangeValue = rangeText
        deleteButton.on_click(self.updatefile)
        cancelButton.on_click(self.cancelDelete)
    
    def cancelDelete(self, b):
        self.delete_menu.children = [self.delete_figures_button]

    def updatefile(self, b):    #Miguel
        
        # load file dsipy
        controller = self.controller
        version = controller.version
        if version == '':
            file_name = controller.name + ".dsipy"
        else:
            file_name = controller.name + '-V' + str(version) + ".dsipy"
        f = open(file_name, "r")
        saved_data = pickle.load(f)
        f.close()

        Figures = saved_data.saved["figure_data"]
        
        #Extract Single Figures
        string = b.indValue.value
        string = string.split(',')
        if b.indValue.value:
            figID = [int(s) for s in string]
        else:
            figID = []
        
        if len(figID) > 0:
            for s in figID:
                del Figures[s-1]
                del controller.figure_data[s-1]
        
        #Extract Range
        string = b.rangeValue.value
        string = string.split('-')
        if b.rangeValue.value:
            indicesRange =[int(s) for s in string]
        else:
            indicesRange = []
        
        if len(indicesRange) == 2:
            del Figures[indicesRange[0]-1:indicesRange[1]]
            del controller.figure_data[indicesRange[0]-1:indicesRange[1]]

        saved_data.saved["figure_data"] = Figures
        
        # controller.figure_data = Figures

        f = open(file_name,"w")
        pickle.dump(saved_data, f)
        f.close()
        sucessMessage = Latex()
        sucessMessage.value = "Figures were deleted and file was updated!"
        self.delete_menu.children = [sucessMessage, self.delete_figures_button]
        
        
        controller.figures.create_figures_widget()
        self.load_widgets()

