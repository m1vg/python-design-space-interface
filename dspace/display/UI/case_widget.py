import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg 


class CaseReport(object):
    
    def __init__(self, controller, by_signature=False):
        setattr(self, 'controller', controller)
        setattr(self, 'by_signature', by_signature)
        
    def create_case_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        by_signature = widgets.CheckboxWidget(
                        description='Cases indicated by signature?', 
                        value=self.by_signature
                        )
        case_id = widgets.TextWidget(description='* Analysis for case:')        
        constraints = widgets.TextareaWidget(description='Biological constraints:')
        button = widgets.ButtonWidget(description='Create Analysis')
        button.on_click(self.create_report)
        button.case_id = case_id
        button.by_signature = by_signature
        button.constraints = constraints
        wi = widgets.ContainerWidget(children=[by_signature,
                                               case_id,
                                               constraints,
                                               button])
        return ('Analyze Case', wi)
        ## controller.update_child('Analyze Case', wi)
    
    def create_report(self, b):
        controller = self.controller
        by_signature = b.by_signature.value
        case_id = str(b.case_id.value)
        constraints = str(b.constraints.value)
        if constraints != '':
            constraints = constraints.split(',')
            constraints = [i.strip() for i in constraints]
        else:
            constraints = None
        display_case = DisplayCase(controller, case_id, by_signature=by_signature, constraints=constraints)
        display_case.create_case_widget()
        
class DisplayCase(object):
    
    def __init__(self, controller, case_id, by_signature=False, pvals=None, constraints=None, subtitle=''):
        setattr(self, 'controller', controller)
        setattr(self, 'pvals', pvals)
        setattr(self, 'case', controller.ds(case_id, 
                                            by_signature=by_signature,
                                            constraints=constraints))
        setattr(self, 'info', None)
        setattr(self, 'equations', None)
        setattr(self, 'log_gains', None)
        setattr(self, 'parameter_table', None)
        setattr(self, 'log_coordinates', False)
        setattr(self, 'dynamic_only', False)
        setattr(self, 'subtitle', subtitle)
        setattr(self, 'fixed_pvals', False if pvals is None else True)
        
    def create_case_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        case = self.case
        self.info = widgets.ContainerWidget()
        self.equations = widgets.ContainerWidget()
        self.log_gains = widgets.ContainerWidget()
        self.parameter_table = widgets.ContainerWidget() 
        self.tolerances_table = widgets.ContainerWidget() 
        save_button = widgets.ButtonWidget(description='Make Nominal Parameter Set')
        save_button.on_click(self.save_parameters)
        if case.is_valid() is True:
            if self.pvals is None:
                self.pvals = case.valid_interior_parameter_set()
        else:
            save_button.visible = False
            self.parameter_table.visible = False
        close_button = widgets.ButtonWidget(description='Close Tab')
        close_button.on_click(self.close_widget)
        wi = widgets.PopupWidget(children=[self.info, 
                                           self.equations,
                                           self.log_gains,
                                           self.parameter_table,
                                           save_button,
                                           self.tolerances_table,
                                           close_button])
        wi.set_css('height', '400px')
        self.update_display()
        subtitle = self.subtitle
        if subtitle != '':
            subtitle = ' (' + subtitle + ')'
        self.title = 'Case ' + self.case.case_number + subtitle
        controller.update_child(self.title, wi)
        
    def close_widget(self, button):
        controller = self.controller 
        controller.update_child(self.title, None)
        
    def save_parameters(self, b):
        controller = self.controller
        controller.pvals.update(self.pvals)

    def update_display(self):
        self.update_info()
        self.update_equations()
        self.update_log_gains()
        self.update_parameter_table()
        self.update_global_tolerances()
        
    def update_info(self):
        controller = self.controller 
        html = widgets.HTMLWidget()
        html_string = '<b>Name:</b> ' + controller.name + '<br>'
        html_string += '<b>Case Number:</b> ' + self.case.case_number + '<br>'
        html_string += '<b>Case Signature:</b> (' + self.case.signature + ')<br>'
        html_string += '<b>Valid:</b> ' + str(self.case.is_valid()) + '<br>'
        html_string += '<b>Is Cyclical:</b> ' + str(self.case.is_cyclical) + '<br><hr>'
        check_box = widgets.CheckboxWidget(description='Logarithmic coordinates', 
                                           value=self.log_coordinates)
        check_box.on_trait_change(self.change_logarithmic, 'value')
        html.value = html_string
        self.info.children = [html, check_box]
        
    def update_equations(self):
        controller = self.controller 
        case = self.case
        check_box = widgets.CheckboxWidget(description='Show only dynamical system?', 
                                           value=self.dynamic_only)
        check_box.on_trait_change(self.change_dynamic_only, 'value')
        if self.case.ssystem.solution is None:
            check_box.disables = True            
        equations = widgets.LatexWidget()
        if self.dynamic_only is False:
            equations.value = case.equations._repr_latex_()
        else:
            ssys = case.ssystem.remove_algebraic_constraints()
            equations.value = ssys.equations._repr_latex_()
        solution = widgets.LatexWidget()
        conditions = widgets.LatexWidget() 
        boundaries = widgets.LatexWidget() 
        if self.log_coordinates is True:
            conditions.value = case.conditions_log._repr_latex_()
            solution_data = case.ssystem.solution_log
            boundary_data = case.boundaries_log
        else:
            conditions.value = case.conditions._repr_latex_()
            solution_data = case.ssystem.solution
            boundary_data = case.boundaries
        if solution_data is not None:
            solution.value = solution_data._repr_latex_()
            boundaries.value = boundary_data._repr_latex_()
        else:
            solution.value = 'S-System has no solution.'
            boundaries.value = 'S-System has no solution.'
        self.equations.children = [check_box,
                                   widgets.HTMLWidget(value='<b>Equations:</b><br>'),
                                   equations, 
                                   widgets.HTMLWidget(value='<b>Conditions:</b><br>'),
                                   conditions, 
                                   widgets.HTMLWidget(value='<b>Solution:</b><br>'),
                                   solution,
                                   widgets.HTMLWidget(value='<b>Boundaries:</b><br>'),
                                   boundaries]

    def update_log_gains(self):
        controller = self.controller 
        case = self.case
        if case.ssystem.solution is None:
            self.log_gains.children = []
            return
        table = widgets.HTMLWidget()
        html_str = '<br><div><table>\n<caption>Logarithmic gains and parameter sensitivities for Case ' + case.case_number + ' (' + case.signature + '). </caption>\n'
        html_str += '<tr ><th rowspan="2" align=center  style="padding:0 15px 0 15px;"> Dependent<br> Variables </th>'
        html_str += '<th colspan="' + str(len(case.independent_variables)) + '" align=center  style="padding:0 15px 0 15px;"> Independent Variables and Parameters</th></tr><tr align=center>'
        for xi in case.independent_variables:
                html_str += '<td><b>{0}</b></td>'.format(xi)
        html_str += '</tr>\n'
        for xd in case.dependent_variables:
            html_str += '<tr><td align=center  style="padding:0 15px 0 15px;"><b>{0}</b></td>'.format(xd)
            for xi in case.independent_variables:
                html_str += '<td align=center  style="padding:0 15px 0 15px;">{0}</td>'.format(str(case.ssystem.log_gain(xd, xi)))
            html_str += '</tr>\n'
        html_str += '</table></div>'
        save_button = widgets.ButtonWidget(description='Save Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        table.value = html_str
        self.log_gains.children = [widgets.HTMLWidget(value='<hr>'),
                                   save_button, 
                                   table,
                                   widgets.HTMLWidget(value='<hr>')]
        return
    
    def update_parameter_table(self):
        controller = self.controller 
        case = self.case
        pvals = self.pvals
        if pvals is None:
            self.parameter_table.children = []
            return
        if case.is_valid(p_bounds=pvals) is False:
            self.parameter_table.children = []
            return            
        table = widgets.HTMLWidget()
        html_str = '<br><div><table>\n<caption>Value for the parameters automatically determined for Case ' + case.case_number + ' (' + case.signature + '). </caption>\n'
        html_str += '<tr ><th align=center  style="padding:0 15px 0 15px;"> Parameters </th><th> Value </th>'
        for xi in sorted(pvals.keys()):
                html_str += '<tr><td><b>{0}</b></td><td>{1}</td></tr>'.format(
                             xi,
                             pvals[xi])
        html_str += '</table></div>'
        table.value = html_str
        save_button = widgets.ButtonWidget(description='Save Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.parameter_table.children = [widgets.HTMLWidget(value='<hr>'),
                                         save_button,
                                         table,
                                         widgets.HTMLWidget(value='<hr>')]
        return
    
    def update_global_tolerances(self):
        controller = self.controller 
        case = self.case
        pvals = self.pvals
        if pvals is None:
            self.tolerances_table.children = []
            return
        if case.is_valid(p_bounds=pvals) is False:
            self.tolerances_table.children = []
            return            
        table = widgets.HTMLWidget()
        html_str = '<br><div><table>\n<caption>Global tolerances determined for Case ' + case.case_number + ' (' + case.signature + ') at the reference parameter values. </caption>\n'
        html_str += '<tr ><th align=center  style="padding:0 15px 0 15px;"> Parameters </th><th> Lower tolerance </th><th> upper tolerance </th>'
        tolerances = case.measure_tolerance(pvals, log_out=self.log_coordinates)
        for xi in sorted(pvals.keys()):
                html_str += '<tr><td><b>{0}</b></td><td>{1}</td><td>{2}</td></tr>'.format(
                             xi,
                             tolerances[xi][0],
                             tolerances[xi][1])
        html_str += '</table></div>'
        table.value = html_str
        save_button = widgets.ButtonWidget(description='Save Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.tolerances_table.children = [widgets.HTMLWidget(value='<hr>'),
                                          save_button,
                                          table,
                                          widgets.HTMLWidget(value='<hr>')]
        return
        
        
    def change_logarithmic(self, name, value):
        self.log_coordinates = value
        self.update_equations()
        
    def change_dynamic_only(self, name, value):
        self.dynamic_only = value
        self.update_equations()
        
    def save_table(self, b):
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
 
       
       
        