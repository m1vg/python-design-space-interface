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
        if 'biological_constraints' not in controller.options:
            bio_constraints = ''
        else:
            bio_constraints = ', '.join(controller.defaults('biological_constraints')) 
        constraints = widgets.TextareaWidget(description='Biological constraints:',
                                             value=bio_constraints)
        button = widgets.ButtonWidget(description='Create Co-Localization')
        button.on_click(self.make_colocalization)
        button.cases = cases
        button.by_signature = by_signature
        button.slice_variables = slice_variables
        button.constraints = constraints
        wi = widgets.ContainerWidget(children=[cases,
                                               by_signature,
                                               slice_variables,
                                               constraints,
                                               button])
        return ('Co-localizations', wi)
        
    def make_colocalization(self, b):
        controller = self.controller
        ds = controller.ds
        case_numbers = str(b.cases.value).split(',')
        case_numbers = [i.strip() for i in case_numbers]
        constraints = str(b.constraints.value)
        if constraints != '':
            constraints = constraints.split(',')
            constraints = [i.strip() for i in constraints if len(i.strip()) > 0]
        else:
            constraints = []
        controller.set_defaults('biological_constraints', constraints)
        cases = ds(case_numbers, by_signature=b.by_signature.value, constraints=constraints)
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
        self.log_coordinates = False
        self.global_tolerance = widgets.ContainerWidget()
        close_button = widgets.ButtonWidget(description='Close Tab')
        close_button.on_click(self.close_widget)
        ss_options = ['log('+ i + ')' for i in controller.ds.dependent_variables]
        dropdown = widgets.DropdownWidget(description='y-axis',
                                          values=ss_options,
                                          value=self.y_variable)
        self.make_plot = widgets.ButtonWidget(description='Create Plot')
        self.make_plot.on_click(self.change_y_axis)
        self.make_plot.yaxis = dropdown
        self.make_plot.visible = False
        dropdown.visible = False
        if len(self.slice_variables) <= 2:
            self.make_plot.visible = True
            if len(self.slice_variables) == 1:
                dropdown.visible = True
        wi = widgets.ContainerWidget(children=[self.info, 
                                               self.constraints_widget,
                                               dropdown,
                                               self.make_plot,
                                               self.global_tolerance,
                                               close_button])
        wi.set_css('height', '400px')
        self.update_display()
        controller.update_child(self.name, wi)
                
    def update_display(self):
        self.update_info()
        self.update_constraints()
        ## self.update_plot()
        self.update_global_tolerances()

        ## self.update_log_gains()
        
    def update_info(self):
        
        title = widgets.HTMLWidget(value='<b> Cases to Co-localize </b>')
        buttons = []
        html_str = '<div><b>Is Valid: {0}</b></div>'.format(self.ci.is_valid())
        if self.ci.is_valid() is False:
            self.make_plot.disabled = True
        else:
            self.make_plot.disabled = False
        valid = widgets.HTMLWidget(value = html_str)
        html_str = '<table><caption> Auxiliary variables for ' + self.name + ' with slice variables ' + ', '.join(self.slice_variables) + '</caption>'
        html_str += '<tr ><th rowspan="2" align=center  style="padding:0 15px 0 15px;"> Slice<br> Variables </th>'
        html_str += '<th colspan="' + str(len(self.cases)) + '" align=center  style="padding:0 15px 0 15px;"> Cases </th></tr><tr align=center>'
        for c in self.cases:
                html_str += '<td style="padding:0 15px 0 15px;"><b>  {0}  </b></td>'.format(c.case_number)
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
        html_str += '</table>'
        html_str += '<caption>Case co-localization assumes that the slice variables '
        html_str += 'for one case are independent from the other cases in the co-localization.'
        html_str += 'Therefore, we define a unique auxiliary variable representing the'
        html_str += 'slice variables for each of cases.</caption>'
        save_table = widgets.ButtonWidget(description='Retain variable table')
        save_table.on_click(self.save_table)
        save_table.table_data = html_str
        variables = widgets.HTMLWidget(value=html_str)
        self.info.children = [title] + buttons + [valid, save_table, variables]
        
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
            ## self.update_plot()
        self.update_constraints()
        
    def update_plot(self):
        controller = self.controller
        ds = controller.ds
        if len(self.slice_variables) > 2:
            return
        ci = self.ci
        if ci.is_valid() is False:
            ## self.plot.children = []
            self.pvals = {}
            return
        pset = self.ci.valid_interior_parameter_set()
        self.pvals = pset[pset.keys()[0]]
        pvals = self.pvals
        fig = plt.figure(dpi=600, facecolor='w')
        cases = [i.case_number for i in self.cases]
        if len(self.slice_variables) == 1:
            fig = plt.figure(figsize=[6, 3], dpi=600, facecolor='w')
            ax1 = fig.add_axes([0.2, 0.3, 0.55, 0.6])
            ax2 = fig.add_axes([0.2, 0.2, 0.55, 0.075])
            ax3 = fig.add_axes([0.8, 0.2, 0.05, 0.7])
            xaxis = self.slice_variables[0]
            xvalues = [pset[i][xaxis] for i in pset]
            x_range = [min(xvalues)*1e-2, max(xvalues)*1e2]
            colors = ds.draw_1D_slice(ax2,
                                      pvals, xaxis, x_range, colorbar=False, intersections=[1, 3, 5, 7])
            ds.draw_1D_ss_function(ax1, self.y_variable,
                                   pvals, xaxis, x_range, colors=colors, lw=2.)
            for case_number in pset:
                case = controller.ds(case_number)
                pvals = pset[case_number]
                yval = case.ssystem.steady_state_function(self.y_variable, pvals)
                ax1.plot(np.log10(pset[case_number][xaxis]), yval, 
                         'o', mfc=colors[case_number], mec='k', ms=5., lw=2.)
            ax1.set_xticklabels('')
            ax1.set_xlabel('')
            ds.draw_region_colorbar(ax3, colors)
            title = 'System design space showing a 1-D case co-localization'
            caption = 'Different cases shown by line segments of different colors. '
            caption += 'Steady state function shown on the y-axis.  '
            caption += 'Bottom bar shows 1D design space slice showing valid cases '
            caption += 'as shown by the colorbar.'
            caption += ' Figure generated with the following parameter values: '
            caption += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals) if i not in [xaxis]]) + '.'
        else:
            options = []
            fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
            ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
            xaxis = self.slice_variables[0]
            yaxis = self.slice_variables[1]
            xvalues = [pset[i][xaxis] for i in pset]
            yvalues = [pset[i][yaxis] for i in pset]
            x_range = [min(xvalues)*1e-2, max(xvalues)*1e2]
            y_range = [min(xvalues)*1e-2, max(xvalues)*1e2]
            colors = ds.draw_2D_slice(ax, pvals, xaxis, yaxis,
                                      x_range, y_range, included_cases=cases)
            for case_number in pset:
                case = controller.ds(case_number)
                pvals = pset[case_number]
                ax.plot(np.log10(pset[case_number][xaxis]), np.log10(pset[case_number][yaxis]), 
                         'o', mfc=colors[case_number], mec='k', ms=5., lw=2.)
            title = 'System design space showing a 2-D case co-localization'
            caption = 'Enumerated co-localized qualitatively-distinct phenotypes represented '
            caption += 'on the z-axis and identified by color.  '
            caption += 'Circles represent automatically determined values for each phenotype.'                
            caption += ' Figure generated with the following parameter values: '
            caption += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals) if i not in [xaxis, yaxis]]) + '.'
        canvas = FigureCanvasAgg(fig) 
        plt.close()
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        ## self.plot.children = options
        controller.figures.add_figure(data, 
                                      title=title,
                                      caption=caption)
        
    
    def update_global_tolerances(self):
        controller = self.controller
        ds = controller.ds
        if len(self.slice_variables) > 2:
            return
        ci = self.ci
        if ci.is_valid() is False:
            self.global_tolerance.children = []
            return
        pvals = self.ci.valid_interior_parameter_set(project=False)
        if pvals is None:
            self.global_tolerance.children = []
            return
        table = widgets.HTMLWidget()
        html_str = '<div><table>\n<caption>Global tolerances determined for ' + self.name + ' showing fold-difference to a large qualitative change{0}. </caption>\n'.format(' in log-coordinates' if self.log_coordinates is True else '') 
        html_str += '<tr ><th align=center  rowspan=2 style="padding:0 15px 0 15px;"> Parameters </th><th colspan=2> Tolerance </th></tr>'
        html_str += '<tr><td style="padding:0 15px 0 15px;"><b> Lower bound</b></td><td style="padding:0 15px 0 15px;"><b> Upper bound</b></td></tr>'
        tolerances = ci.measure_tolerance(pvals, log_out=self.log_coordinates)
        for xi in sorted(pvals.keys()):
            lower_th = 1e-15 if self.log_coordinates is False else -15
            upper_th = 1e15 if self.log_coordinates is False else 15
            lower, upper = tolerances[xi]
            html_str += '<tr><td style="padding:0 15px 0 15px;"><b>{0}</b></td><td style="padding:0 15px 0 15px;">{1}</td><td style="padding:0 15px 0 15px;">{2}</td></tr>'.format(
                         xi,
                         lower if lower > lower_th else '-&infin;',
                         upper if upper < upper_th else '&infin;')
        html_str += '</table><caption>'
        html_str += 'Note: Global tolerance calculated based on the following values for the parameters: ' 
        html_str += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals.keys())]) + '.'
        html_str += '</caption></div>'
        table.value = html_str
        save_button = widgets.ButtonWidget(description='Retain Global Tolerance Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.global_tolerance.children = [widgets.HTMLWidget(value='<br>'),
                                          save_button,
                                          table]
        return
    
    def save_table(self, b):
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)

    def change_y_axis(self, b):
        self.y_variable = str(b.yaxis.value)
        ## self.plot.children = [self.plot.children[0]]
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

        
        
