import dspace
import dspace.plotutils
import dspace.display

import numpy as np
from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib as mt
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from subprocess import call, Popen, PIPE
from dspace.graphs.designspace_graph import GraphGenerator
from dspace.display.UI.case_widget import DisplayCase
import base64

from collections import OrderedDict

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
        self.log_coordinates = True
        self.global_tolerance = widgets.ContainerWidget()
        close_button = widgets.ButtonWidget(description='Close Tab')
        close_button.on_click(self.close_widget)
        ss_options = ['log('+ i + ')' for i in controller.ds.dependent_variables]
        self.y_dropdown = widgets.DropdownWidget(description='y-axis',
                                          values=ss_options,
                                          value=self.y_variable)
        self.make_plot = widgets.ButtonWidget(description='Create Plot')
        self.make_plot.on_click(self.change_y_axis)
        self.make_plot.yaxis = self.y_dropdown
        self.make_plot.visible = True
        check_box = widgets.CheckboxWidget(description='Logarithmic coordinates', 
                                           value=self.log_coordinates)
        check_box.on_trait_change(self.update_log, 'value')
        self.y_dropdown.visible = False
        if len(self.slice_variables) <= 3:
            ## self.make_plot.visible = True
            if len(self.slice_variables) == 1:
                self.y_dropdown.visible = True
            if len(self.slice_variables) == 3:
                self.y_dropdown = widgets.HTMLWidget(value='<font color="red">Warning 3D plots are experimental.</font>')
        wi = widgets.ContainerWidget(children=[self.info, 
                                               self.constraints_widget,
                                               check_box,
                                               self.y_dropdown,
                                               self.make_plot,
                                               self.global_tolerance,
                                               close_button])
        wi.set_css('height', '400px')
        self.update_display()
        controller.update_child(self.name, wi)
                
    def update_display(self):
        self.update_info()
        self.update_constraints()
        self.update_global_tolerances()
        
    def update_log(self, name, value):
        controller = self.controller
        self.log_coordinates = value
        if value == False:
            ss_old = ['log('+ i + ')' for i in controller.ds.dependent_variables]
            ss_new = [i for i in controller.ds.dependent_variables]
        else:
            ss_old = [i for i in controller.ds.dependent_variables]
            ss_new = ['log('+i+')' for i in controller.ds.dependent_variables]
        ss_options = OrderedDict([(unicode(i),i) for i in ss_new])
        index = ss_old.index(self.y_variable)
        self.y_variable = ss_new[index]
        self.y_dropdown.values = ss_options
        self.y_dropdown.value = unicode(self.y_variable)
        self.update_display()
        
        
    def update_info(self):
        
        title = widgets.HTMLWidget(value='<b> Cases to Co-localize </b>')
        buttons = []
        html_str = '<div><b>Is Valid: {0}</b></div>'.format(self.ci.is_valid())
        if self.ci.is_valid() is False:
            self.make_plot.disabled = True
        else:
            self.make_plot.disabled = False
        valid = widgets.HTMLWidget(value = html_str)
        html_str = '<table><caption> Auxiliary variables for ' + self.name 
        html_str += ' with ' + ', '.join(self.slice_variables) 
        html_str += ' as the slice variable{0}.</caption>'.format('s' if len(self.slice_variables) > 1 else '')
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
        html_str += 'for one case are independent from the slice variables for the other cases in the co-localization.'
        html_str += 'Each auxiliary variable corresponds to a slice variable for one cases.</caption>'
        save_table = widgets.ButtonWidget(description='Save variable table')
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
        ## if len(self.slice_variables) > 3:
        ##     return
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
        elif len(self.slice_variables) == 2:
            options = []
            fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
            ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
            xaxis = self.slice_variables[0]
            yaxis = self.slice_variables[1]
            xvalues = [pset[i][xaxis] for i in pset]
            yvalues = [pset[i][yaxis] for i in pset]
            x_range = [min(xvalues)*1e-2, max(xvalues)*1e2]
            y_range = [min(yvalues)*1e-2, max(yvalues)*1e2]
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
        elif len(self.slice_variables) == 3:
            options = []
            fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
            ax = fig.add_axes([0.2, 0.2, 0.7, 0.7], projection='3d')
            xaxis = self.slice_variables[0]
            yaxis = self.slice_variables[1]
            zaxis = self.slice_variables[2]
            xvalues = [pset[i][xaxis] for i in pset]
            yvalues = [pset[i][yaxis] for i in pset]
            zvalues = [pset[i][zaxis] for i in pset]
            x_range = [min(xvalues)*1e-2, max(xvalues)*1e2]
            y_range = [min(yvalues)*1e-2, max(yvalues)*1e2]
            z_range = [min(zvalues)*1e-2, max(zvalues)*1e2]
            colors = ds.draw_3D_slice(ax, pvals, xaxis, yaxis, zaxis,
                                      x_range, y_range, z_range, included_cases=cases)
            title = 'System design space showing a 3-D case co-localization'
            caption = 'Enumerated co-localized qualitatively-distinct phenotypes represented '
            caption += 'on the z-axis and identified by color.  '
            caption += ' Figure generated with the following parameter values: '
            caption += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals) if i not in [xaxis, yaxis]]) + '.'
        else:
            options = []
            slice_variables = self.slice_variables
            fig = plt.figure(figsize=[6, 4], dpi=600, facecolor='w')
            ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
            values = {i:[pset[i][axis] for axis in slice_variables] for i in pset}
            ## ranges = {axis:(min(values[axis])*1e-2,max(values[axis])*1e2) for axis in slice_variables}
            colors = {name:mt.cm.hsv(float(i)/float(len(cases))) for i, name in enumerate(cases)}
            min_value = None
            max_value = None
            for case in cases: 
                y_values = np.log10(values[case])
                min_y = min(y_values)
                max_y = min(y_values)
                if min_value is None:
                    min_value = min_y
                else:
                    min_value = min([max_value, max_y])
                if max_value is None:
                    max_value = min_y
                else:
                    max_value = max([max_value, max_y])
                ax.plot(range(len(slice_variables)), 
                        y_values,
                        lw=2., c=colors[case])
            ax.set_ylim([min_value, max_value])
            ax.set_xlim([0, len(slice_variables)])
            ax.set_xticks(range(len(slice_variables)))
            ax.set_xticklabels(['$' + controller.symbols[i] + '$' for i in slice_variables])
            title = 'Values for the slice variable for the n-D case co-localization'
            caption = 'The y-axis represents value for the slice variable on the'
            caption += ' x-axis for a case identified by color.'
            caption += ' Figure generated with the following parameter values: '
            caption += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals) if i[0] != '$']) + '.'
        canvas = FigureCanvasAgg(fig) 
        plt.close()
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        controller.figures.add_figure(data, 
                                      title=title,
                                      caption=caption)
        
    
    def update_global_tolerances(self):
        controller = self.controller
        ds = controller.ds
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
        save_button = widgets.ButtonWidget(description='Save Global Tolerance Table')
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
        controller.save_widget_data(b)

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

        
        
