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
from dspace.display.UI.case_widget import DisplayCase
import base64

class CaseColocalization(object):
    
    def __init__(self, controller, by_signature=False):
        setattr(self, 'controller', controller)
        setattr(self, 'colocalization_data', widgets.ContainerWidget())
        setattr(self, 'by_signature', by_signature)
    
    def colocalization_widget(self):
        controller = self.controller
        cases = widgets.TextareaWidget(description='* Cases to co-localize:')
        by_signature = widgets.CheckboxWidget(description='Cases indicated by signature?', 
                                              value=self.by_signature)
        slice_variables = widgets.TextareaWidget(description='* Slice variables:')
        button = widgets.ButtonWidget(description='Create Co-Localization')
        button.on_click(self.make_colocalization)
        button.cases = cases
        button.by_signature = by_signature
        button.slice_variables = slice_variables
        wi = widgets.ContainerWidget(children=[cases,
                                               by_signature,
                                               slice_variables,
                                               button])
        return ('Co-localizations', wi)
        
    def make_colocalization(self, b):
        controller = self.controller
        ds = controller.ds
        case_numbers = str(b.cases.value).split(',')
        case_numbers = [i.strip() for i in case_numbers]
        cases = ds(case_numbers, by_signature=b.by_signature.value)
        slice_variables = str(b.slice_variables.value).split(',')
        slice_variables = [i.strip() for i in slice_variables]
        co = DisplayColocalization(self.controller, cases, slice_variables)
        co.make_display_widget()
        
        
class DisplayColocalization(object):
    
    def __init__(self, controller, cases, slice_variables):
        setattr(self, 'controller', controller)
        setattr(self, 'colocalizations_widget', None)
        setattr(self, 'cases', cases)
        setattr(self, 'slice_variables', slice_variables)
        setattr(self, 'name', 'Co-localization: ' + ', '.join([i.case_number for i in self.cases]))
        setattr(self, 'constraints', [])
        setattr(self, 'active_constraints', False)
        setattr(self, 'ci', dspace.CaseColocalization(cases, slice_variables))
        setattr(self, 'pvals', {})
        setattr(self, 'y_variable', 'log('+controller.ds.dependent_variables[0]+')')
        
    def make_display_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        self.info = widgets.ContainerWidget()
        self.constraints_widget = widgets.ContainerWidget()
        self.plot = widgets.ContainerWidget()
        close_button = widgets.ButtonWidget(description='Close Tab')
        close_button.on_click(self.close_widget)
        wi = widgets.ContainerWidget(children=[self.info, 
                                               self.constraints_widget,
                                               self.plot,
                                               close_button])
        wi.set_css('height', '400px')
        self.update_display()
        controller.update_child(self.name, wi)
                
    def update_display(self):
        self.update_info()
        self.update_constraints()
        self.update_plot()

        ## self.update_log_gains()
        
    def update_info(self):
        
        title = widgets.HTMLWidget(value='<b> Cases to Co-localize </b>')
        buttons = []
        html_str = '<div><b>Is Valid: {0}</b></div>'.format(self.ci.is_valid())
        html_str += '<br><div><table><caption> Auxiliary variables for co-localized cases.'
        html_str += '<tr ><th rowspan="2" align=center  style="padding:0 15px 0 15px;"> Slice<br> Variables </th>'
        html_str += '<th colspan="' + str(len(self.cases)) + '" align=center  style="padding:0 15px 0 15px;"> Cases </th></tr><tr align=center>'
        for c in self.cases:
                html_str += '<td><b>  {0}  </b></td>'.format(c.case_number)
        html_str += '</tr>\n'
        pset = self.ci.valid_interior_parameter_set()
        for i in self.cases:
            key = i.case_number
            case_button = widgets.ButtonWidget(description='Case ' + key)
            buttons.append(case_button)
            case_button.pvals = pset[key] if key in pset else None
            case_button.on_click(self.open_case)
        for j, xj in enumerate(self.slice_variables):
            html_str += '<tr align=center><td>{0}</td>'.format(xj)
            for i, c in enumerate(self.cases):
                html_str += '<td>${0}_{1}</td>'.format(xj, i)
            html_str += '</tr>'                
        variables = widgets.HTMLWidget(value=html_str)
        self.info.children = [title] + buttons + [variables]
        
    def update_constraints(self):
        constraints_widget = widgets.TextareaWidget(description='Constraints',
                                                    value = ',\n'.join(self.constraints)
                                                    )
        constraints_widget.visible = self.active_constraints
        button = widgets.ButtonWidget(description='Done' if self.active_constraints else 'Modify constraints')
        button.constraints_widget = constraints_widget
        button.on_click(self.modify_constraints)
        self.constraints_widget.children = [constraints_widget, button]
        
    def modify_constraints(self, b):
        if self.active_constraints is False:
            self.active_constraints = True
        else:
            self.active_constraints = False
            self.constraints = str(b.constraints_widget.value).split(',')
            self.constraints = [i.strip() for i in self.constraints]
            self.ci = dspace.CaseColocalization(self.cases, 
                                                self.slice_variables,
                                                constraints = self.constraints)
            self.update_info()
            self.update_plot()
        self.update_constraints()
        
    def update_plot(self):
        controller = self.controller
        ds = controller.ds
        if len(self.slice_variables) > 2:
            return
        ci = self.ci
        if ci.is_valid() is False:
            self.plot.children = []
            self.pvals = {}
            return
        pset = self.ci.valid_interior_parameter_set()
        self.pvals = pset[pset.keys()[0]]
        pvals = self.pvals
        fig = plt.figure(dpi=600, facecolor='w')
        cases = [i.case_number for i in self.cases]
        if len(self.slice_variables) == 1:
            ss_options = ['log('+ i + ')' for i in controller.ds.dependent_variables]
            dropdown = widgets.DropdownWidget(description='y-axis', 
                                              values=ss_options,
                                              value=self.y_variable)
            dropdown.on_trait_change(self.change_y_axis, 'value')
            options = [dropdown]
            ax1 = fig.add_axes([0.1, 0.2, 0.7, 0.7])
            ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.05])
            ax3 = fig.add_axes([0.85, 0.1, 0.05, 0.8])
            xaxis = self.slice_variables[0]
            xvalues = [pset[i][xaxis] for i in pset]
            x_range = [min(xvalues)*1e-1, max(xvalues)*1e1]
            colors = ds.draw_1D_slice(ax2,
                                      pvals, xaxis, x_range, colorbar=False)
            ds.draw_1D_ss_function(ax1, self.y_variable,
                                   pvals, xaxis, x_range, colors=colors)
            ax1.set_xticklabels('')
            ds.draw_region_colorbar(ax3, colors)

        else:
            options = []
            ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
            xaxis = self.slice_variables[0]
            yaxis = self.slice_variables[1]
            xvalues = [pset[i][xaxis] for i in pset]
            yvalues = [pset[i][yaxis] for i in pset]
            x_range = [min(xvalues)*1e-1, max(xvalues)*1e1]
            y_range = [min(xvalues)*1e-1, max(xvalues)*1e1]
            ds.draw_2D_slice(ax, pvals, xaxis, yaxis, x_range, y_range, included_cases=cases)            
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        image_widget = widgets.ImageWidget(value=data)
        image_widget.set_css('height', '400px')
        image_widget.set_css('width', '480px')
        self.plot.children = options + [image_widget]
        plt.close()
        
    def change_y_axis(self, name, value):
        self.y_variable = str(value)
        self.plot.children = [self.plot.children[0]]
        self.update_plot()
        
    def open_case(self, b):
        case_id = b.description.split(' ')[1]
        c_widget = DisplayCase(self.controller, 
                               str(case_id), 
                               by_signature=False,
                               pvals=b.pvals,
                               subtitle='Co-localized')
        c_widget.create_case_widget()

        
    def close_widget(self, b):
        controller = self.controller 
        controller.update_child(self.name, None)

        
        
