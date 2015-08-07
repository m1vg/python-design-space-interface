import dspace
import dspace.plotutils
import dspace.display

import numpy as np
from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib as mt
import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from subprocess import call, Popen, PIPE
from dspace.graphs.designspace_graph import GraphGenerator

import base64

def eigenvalue_compare(eigenvalues, component='real', rank=1):
    if component == 'real':
        eig = eigenvalues.real
    else:
        eig = eigenvalues.imag
    value = sorted(eig)
    rank = min(rank, len(eig))
    return value[-rank]

class MakePlot(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'plot_data', widgets.ContainerWidget())
        setattr(self, 'title', None)
        setattr(self, 'caption', None)
    
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
        cmd = Popen(['dot'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = cmd.communicate(input='')        
        if len(err) == 0:
            widget_types.append('Network Graph')
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
        xlabel = widgets.DropdownWidget(
                  description='* X-Axis',
                  values=controller.ds.independent_variables, 
                  value=xaxis)
        ylabel = widgets.DropdownWidget(
                  description='* Y-Axis',
                  values=controller.ds.independent_variables,
                  value=yaxis)
        xmin = widgets.FloatTextWidget(description='* X-Min',
                                       value=range_x[0])
        xmax = widgets.FloatTextWidget(description='* X-Max',
                                       value=range_x[1])
        ymin = widgets.FloatTextWidget(description='* Y-Min',
                                       value=range_y[0])
        ymax = widgets.FloatTextWidget(description='* Y-Max',
                                       value=range_y[1])
        center_axes = widgets.CheckboxWidget(description='Center Axes',
                                             value=center)
        plot_type = widgets.DropdownWidget(description='* Plot Type',
                                           values=self.widget_types,
                                           value='Design Space (interactive)')
        title_widget = widgets.TextWidget(description='Title')
        caption_widget = widgets.TextareaWidget(description='Caption')
        included = controller.defaults('included_cases')
        if included is None:
            included = []
        if isinstance(included, list) is False:
            included = [included]
        included_widget = widgets.TextareaWidget(description='Only Cases',
                                                 value=', '.join(included))
        wi = widgets.ContainerWidget(children=[xlabel, ylabel, plot_type, 
                                               xmin, xmax, ymin, ymax,
                                               center_axes, title_widget, 
                                               caption_widget, included_widget])
        plot_type.widget_container = wi
        plot_type.on_trait_change(self.update_plot_widget, 'value')
        button = widgets.ButtonWidget(value=False, description='Add Plot')
        button.on_click(self.make_plot)
        button.xlabel = xlabel
        button.ylabel = ylabel
        button.xmin = xmin
        button.xmax = xmax
        button.ymin = ymin
        button.ymax = ymax
        button.center_axes = center_axes
        button.plot_type = plot_type
        button.title = title_widget
        button.caption = caption_widget
        button.included = included_widget
        button.wi = wi
        self.title = title_widget
        self.caption = caption_widget
        add_plot = widgets.ContainerWidget(description='Add Plot', children=[wi, self.plot_data, button])
        controller.update_child('Create Plot', add_plot)
        self.update_plot_widget('value', 'Design Space (Interactive)')
        
    def update_plot_widget(self, name, value):
        controller = self.controller
        zlim = controller.defaults('zlim')
        if value == 'Design Space (interactive)':
            wi = widgets.ContainerWidget(children=[])
            self.plot_data.children = [wi]
            self.title.value = 'System design space'
            self.caption.value = 'System design space with the enumerated qualitatively-distinct phenotypes represented on the z-axis, identified by color.'                
        elif value == 'Design Space':
            intersections_widget = widgets.DropdownWidget(description='# Intersetcions', 
                                                          values=['Single',
                                                                  'Single and Three',
                                                                  'All',],
                                                          value='Single and Three')
            wi = widgets.ContainerWidget(children=[intersections_widget])
            wi.intersections = intersections_widget
            self.title.value = 'System design space'
            self.caption.value = 'Enumerated qualitatively-distinct phenotypes represented on the z-axis and identified by color.'                
            self.plot_data.children = [wi]
        elif value == 'Stability':
            resolution_widget = widgets.FloatTextWidget(description='Resolution', value=100)
            wi = widgets.ContainerWidget(children=[resolution_widget])
            wi.resolution = resolution_widget
            self.plot_data.children = [wi]
            self.title.value = 'System design space showing stability of the fixed points'
            self.caption.value = 'Number of eigenvalues with positive real part represented as a heat map on the z-axis.'
        elif value == 'Eigenvalues':
            component_widget = widgets.DropdownWidget(description='Complex component',
                                                      values=['Real', 'Imaginary'],
                                                      value='Real')
            resolution_widget = widgets.FloatTextWidget(description='Resolution', value=100)
            zlim_auto = (zlim is None)
            zlim_widget = widgets.CheckboxWidget(description='Automatic Z-Lim', value=zlim_auto)
            if zlim_auto is True:
                zlim = [0., 0.]
            zmin_widget = widgets.FloatTextWidget(description='Z-Min', value=zlim[0])
            zmax_widget = widgets.FloatTextWidget(description='Z-Max', value=zlim[1])
            parallel_widget = widgets.CheckboxWidget(description='Compute in Parallel', value=False)
            number_dynamic = len(controller.ds.dependent_variables)
            number_dynamic -= len(controller.ds.auxiliary_variables)
            select_widget = widgets.DropdownWidget(description='Rank to Plot',
                                                   values = [str(i+1) for i in range(number_dynamic)],
                                                   value=str(1))
            wi = widgets.ContainerWidget(children=[component_widget, 
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
        elif value in ['Steady State Concentration', 'Steady State Flux', 'Steady State Function']:
            log_linear_widget = widgets.CheckboxWidget(description='Function is log linear',
                                                       value=True)
            if value == 'Steady State Flux':
                flux_options = ['log(V_'+ i + ')' for i in controller.ds.dependent_variables]
                function_widget = widgets.DropdownWidget(values=flux_options)
                self.title.value = 'System design space showing a steady state flux'
                self.caption.value = 'Steady state flux shown as a heat map on the z-axis.'
            elif value == 'Steady State Function':
                function_widget = widgets.TextWidget(description='Function', value='')
                log_linear_widget.value = False
                self.title.value = 'System design space showing a function at steady state'
                self.caption.value = 'Steady state function shown as a heat map on the z-axis.'
            else:
                ss_options = ['log('+ i + ')' for i in controller.ds.dependent_variables]
                function_widget = widgets.DropdownWidget(values=ss_options)
                self.title.value = 'System Design Space showing a steady state concentration'
                self.caption.value = 'Steady state concentration shown as a heat map on the z-axis.'
            resolution_widget = widgets.FloatTextWidget(description='Resolution', value=100)
            parallel_widget = widgets.CheckboxWidget(description='Compute in Parallel', value=False)
            zlim_auto = (zlim is None)
            zlim_widget = widgets.CheckboxWidget(description='Automatic Z-Lim', value=zlim_auto)
            if zlim_auto is True:
                zlim = [0., 0.]
            zmin_widget = widgets.FloatTextWidget(description='Z-Min', value=zlim[0])
            zmax_widget = widgets.FloatTextWidget(description='Z-Max', value=zlim[1])
            wi = widgets.ContainerWidget(children=[function_widget, resolution_widget,
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
        elif value =='Network Graph':
            type_widget = widgets.DropdownWidget(description='Layout Type',
                                                 values=['Hierarchical',
                                                         'Radial',
                                                         'Circular'],
                                                 value='Hierarchical')
            variables_widget = widgets.TextareaWidget(description='Include Variables',
                                                      value='')
            using_pvals = widgets.CheckboxWidget(description='Using Parameters',
                                                 value=False)
            show_r = widgets.CheckboxWidget(description='Show Regulation',
                                            value=True)
            wi = widgets.ContainerWidget(children=[type_widget, variables_widget,
                                                   using_pvals, show_r])
            wi.type_w = type_widget
            wi.variables = variables_widget
            wi.using_pvals = using_pvals
            wi.show = show_r
            self.plot_data.children = [wi]
            self.title.value = 'Automatically generated network architecture showing mass flow and regulatory interactions'
            self.caption.value = 'Network of interactions within the system. Mass flow between species (black arrows); information flow affecting processes (gray arrows).'
        if controller.name != '':
            title = 'Analysis of the ' + controller.name + ' by ' + self.title.value.lower()
            self.title.value = title
            
    def make_plot(self, b):
        controller = self.controller
        b.description = 'Creating plot... Please Wait.'
        b.disabled = True
        if b.plot_type.value == 'Design Space (interactive)':
            self.make_interactive_plot(b)
        elif b.plot_type.value == 'Design Space':
            self.make_static_plot(b)
        elif b.plot_type.value == 'Stability':
            self.make_stability_plot(b)
        elif b.plot_type.value == 'Eigenvalues':
            self.make_eigenvalue_plot(b)
        elif b.plot_type.value == 'Network Graph':
            self.make_network_graph(b)
        else:
            self.make_function_plot(b)
        b.description = 'Add Plot'
        b.disabled = False
    
    def axes_ranges(self, b):
        pvals = self.controller.pvals
        ranges = [[b.xmin.value, b.xmax.value],[b.ymin.value, b.ymax.value]]
        if b.center_axes.value is False:
            return ranges
        ranges = [[pvals[str(b.xlabel.value)]*i for i in ranges[0]],
                  [pvals[str(b.ylabel.value)]*i for i in ranges[1]]]
        return ranges
    
    def included_cases(self, b):
        included = str(b.included.value).strip()
        if len(included) == 0:
            return None
        included = [i.strip() for i in included.split(',')]
        return included
        
    def make_interactive_plot(self, b):
        controller = self.controller
        button = widgets.ButtonWidget(description='Stop interactive plot')
        button.on_click(self.remove_plot)
        button.name = 'Interactive Plot (' + str(np.random.randint(0, 1000)) + ')'
        image_widget = widgets.ImageWidget()
        rangex, rangey = self.axes_ranges(b)
        interactive_plot = controller.ds.draw_2D_slice_notebook(controller.pvals, str(b.xlabel.value),
                                                                str(b.ylabel.value),
                                                                rangex, rangey,
                                                                {i:[1e-5, 1e5] for i in
                                                                 controller.pvals},
                                                                intersections=[1],
                                                                image_container=image_widget)
        wi = widgets.ContainerWidget(description=button.name,
                                     children=[interactive_plot, 
                                               button,
                                               image_widget])
        controller.update_child(button.name, wi)
        
    def make_static_plot(self, b):
        controller = self.controller
        fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        ax.set_title('Design Space plot')
        plot_data = self.plot_data.children[0]
        intersects = plot_data.intersections.value
        intersections_dict = {'Single':[1],
                              'Single and Three':[1,3],
                              'All':range(1, 100)}
        rangex, rangey = self.axes_ranges(b)
        controller.ds.draw_2D_slice(ax, controller.pvals, str(b.xlabel.value), str(b.ylabel.value),
                                    rangex, rangey,
                                    intersections=intersections_dict[intersects],
                                    included_cases=self.included_cases(b))
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value)
        fig=plt.clf()
        
    def make_stability_plot(self, b):
        controller = self.controller
        fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        ax.set_title('Stability plot')
        plot_data = self.plot_data.children[0]
        resolution = plot_data.resolution.value
        rangex, rangey = self.axes_ranges(b)
        controller.ds.draw_2D_positive_roots(ax, controller.pvals, str(b.xlabel.value), str(b.ylabel.value),
                                             rangex, rangey,
                                             resolution=resolution,
                                             included_cases=self.included_cases(b))
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value)
        fig=plt.clf()
    
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
        controller.ds.draw_2D_ss_function(ax, function, controller.pvals, 
                                          str(b.xlabel.value),
                                          str(b.ylabel.value),
                                          rangex, rangey, zlim=zlim,
                                          log_linear=log_linear, resolution=resolution, 
                                          parallel=parallel,
                                          included_cases=self.included_cases(b))
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value)
        fig=plt.clf()
        
    def make_eigenvalue_plot(self, b):
        controller = self.controller
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
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value)
        fig=plt.clf()
        
    def make_network_graph(self, b):
        controller = self.controller
        plot_data = self.plot_data.children[0]
        graph_type = {'Hierarchical':'dot',
                      'Radial':'twopi',
                      'Circular':'circo'}
        g=GraphGenerator(controller.ds)
        variables = str(plot_data.variables.value)
        variables = [i.strip() for i in variables.split(',')]
        variables = [i for i in variables if i in controller.ds.independent_variables]
        use_pvals = plot_data.using_pvals.value
        if use_pvals is True:
            pvals = controller.pvals
        else:
            pvals = None
        fig = plt.figure(dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        rangex, rangey = self.axes_ranges(b)
        ax.set_title('Network Graph')
        controller.ds.draw_network_graph(ax,
                                         p_vals=pvals,
                                         graph_type=graph_type[str(plot_data.type_w.value)],
                                         included_variables=variables,
                                         cmap=mt.cm.hot_r,
                                         show_regulation=plot_data.show.value)
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, title=b.title.value, caption=b.caption.value)
        fig=plt.clf()
        
    def remove_plot(self, b):
        controller = self.controller
        controller.update_child(b.name, None)
        
        
class DisplayFigures(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'figures_widget', None)
        
    def create_figures_widget(self):
        
        controller = self.controller
        self.figures_widget = widgets.ContainerWidget()
        controller.update_child('Figures', self.figures_widget)
        
    def add_figure(self, image_data, title='', caption = ''):
        controller = self.controller
        figures = controller.figure_data
        figures.append((image_data, title, caption))
        self.add_figure_widget(image_data, title=title, caption = caption)
        
    def add_figure_widget(self, image_data, title='', caption = ''):
        image_widget = widgets.ImageWidget()
        image_widget.value = image_data
        children = [i for i in self.figures_widget.children]      
        if len(title) > 0:
            title = title + '.'
        if len(caption) > 0:
            caption = '  ' + caption
        html_str = '<b>Figure '+str(len(children)+1)+'.  '+title+'</b>' + caption
        html_widget = widgets.HTMLWidget(value=html_str)
        wi = widgets.PopupWidget(children=[image_widget, html_widget])
        wi.set_css('height', '360px')
        wi.set_css('width', '360px')
        children.append(wi)
        self.figures_widget.children = children
        
    def load_widgets(self):
        controller = self.controller
        for data in controller.figure_data:
            self.add_figure_widget(data[0], title=data[1], caption=data[2])
        
        