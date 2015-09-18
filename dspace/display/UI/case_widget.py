import dspace
import dspace.plotutils
import dspace.display

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
    VBox = Box
    HBox = Box
    ipy_old = True
else:
    from ipywidgets import *
    from popup import Popup as PopupWidget
    def Popup(children=[], **kwargs):
        return PopupWidget(children=[VBox(children=children)], **kwargs)
    ipy_old = False
    
from IPython.display import clear_output, display

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
        by_signature = Checkbox(
                        description='Cases indicated by signature?', 
                        value=self.by_signature
                        )
        case_id = Text(description='* Analysis for case:')        
        if 'biological_constraints' not in controller.options:
            bio_constraints = ''
        else:
            bio_constraints = ', '.join(controller.defaults('biological_constraints')) 
        constraints = Textarea(description='Biological constraints:',
                                             value=bio_constraints)
        button = Button(description='Create Analysis')
        button.on_click(self.create_report)
        button.case_id = case_id
        button.by_signature = by_signature
        button.constraints = constraints
        wi = VBox(children=[by_signature,
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
            constraints = [i.strip() for i in constraints if len(i.strip()) > 0]
        else:
            constraints = []
        controller.set_defaults('biological_constraints', constraints)
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
        ## setattr(self, 'fixed_pvals', False if pvals is None else True)
        
    def create_case_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        case = self.case
        self.info = VBox()
        self.equations = VBox()
        self.log_gains = VBox()
        self.parameter_table = VBox() 
        self.tolerances_table = VBox() 
        self.bounding_box_table = VBox()
        calculate_pvals = Button(description='Determine values for the parameters')
        calculate_pvals.visible = False
        calculate_pvals.on_click(self.identify_parameters)
        if case.is_valid() is True:
            if self.pvals is None:
                calculate_pvals.visible = True
        close_button = Button(description='Close Tab')
        close_button.on_click(self.close_widget)
        wi = Popup(children=[self.info, 
                                           self.equations,
                                           self.log_gains,
                                           self.bounding_box_table,
                                           calculate_pvals,
                                           self.parameter_table,
                                           self.tolerances_table,
                                           close_button])
        if ipy_old is True:
            wi.set_css('height', '400px')
        else:
            wi.width = '700px'
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
        self.update_bounding_box()
        
    def update_info(self):
        controller = self.controller 
        html = HTML()
        html_string = '<b>Name:</b> ' + controller.name + '<br>'
        html_string += '<b>Case Number:</b> ' + self.case.case_number + '<br>'
        html_string += '<b>Case Signature:</b> (' + self.case.signature + ')<br>'
        html_string += '<b>Valid:</b> ' + str(self.case.is_valid()) + '<br>'
        html_string += '<b>Is Cyclical:</b> ' + str(self.case.is_cyclical) + '<br><hr>'
        check_box = Checkbox(description='Logarithmic coordinates', 
                                           value=self.log_coordinates)
        check_box.on_trait_change(self.change_logarithmic, 'value')
        html.value = html_string
        self.info.children = [html, check_box]
        
    def update_equations(self):
        controller = self.controller 
        case = self.case
        check_box = Checkbox(description='Show only dynamical system?', 
                                           value=self.dynamic_only)
        check_box.on_trait_change(self.change_dynamic_only, 'value')
        if self.case.ssystem.solution is None:
            check_box.disables = True            
        equations = Latex()
        if self.dynamic_only is False:
            equations.value = case.equations._repr_latex_()
        else:
            ssys = case.ssystem.remove_algebraic_constraints()
            equations.value = ssys.equations._repr_latex_()
        solution = Latex()
        conditions = Latex() 
        boundaries = Latex() 
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
                                   HTML(value='<b>Equations:</b><br>'),
                                   equations, 
                                   HTML(value='<b>Conditions:</b><br>'),
                                   conditions, 
                                   HTML(value='<b>Solution:</b><br>'),
                                   solution,
                                   HTML(value='<b>Boundaries:</b><br>'),
                                   boundaries]

    def update_log_gains(self):
        controller = self.controller 
        case = self.case
        if case.ssystem.solution is None:
            self.log_gains.children = []
            return
        table = HTML()
        html_str = '<div><table>\n<caption>Logarithmic gains and parameter sensitivities for Case ' + case.case_number + ' (' + case.signature + '). </caption>\n'
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
        html_str += '</table></div><br>'
        save_button = Button(description='Save Log-Gain Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        table.value = html_str
        self.log_gains.children = [HTML(value='<br>'),
                                   save_button, 
                                   table]
        return
    
    def update_bounding_box(self):
        controller = self.controller 
        case = self.case
        pvals = self.pvals
        if case.is_valid(p_bounds=pvals) is False:
            self.bounding_box_table.children = []
            return            
        table = HTML()
        html_str = '<div><table>\n<caption>Bounding box determined for Case ' + case.case_number + ' (' + case.signature + '){0}. </caption>\n'.format(' in log-coordinates' if self.log_coordinates is True else '') 
        html_str += '<tr ><th align=center  rowspan=2 style="padding:0 15px 0 15px;"> Parameters </th><th colspan=2> Tolerance </th></tr>'
        html_str += '<tr><td style="padding:0 15px 0 15px;"><b> Lower bound</b></td><td style="padding:0 15px 0 15px;"><b> Upper bound</b></td></tr>'
        tolerances = case.bounding_box(log_out=self.log_coordinates)
        for xi in sorted(tolerances.keys()):
            lower_th = 1e-15 if self.log_coordinates is False else -15
            upper_th = 1e15 if self.log_coordinates is False else 15
            lower, upper = tolerances[xi]
            html_str += '<tr><td style="padding:0 15px 0 15px;"><b>{0}</b></td><td style="padding:0 15px 0 15px;">{1}</td><td style="padding:0 15px 0 15px;">{2}</td></tr>'.format(
                         xi,
                         lower if lower > lower_th else '-&infin;',
                         upper if upper < upper_th else '&infin;')
        html_str += '</table></div>'
        table.value = html_str
        save_button = Button(description='Save Bounding Box Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.bounding_box_table.children = [HTML(value='<br>'),
                                            save_button,
                                            table]
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
        make_nominal = Button(description='Make Nominal Parameter Set')
        make_nominal.on_click(self.save_parameters)
        table = HTML()
        html_str = '<div><table>\n<caption>Value for the parameters automatically determined for Case ' + case.case_number + ' (' + case.signature + '). </caption>\n'
        html_str += '<tr ><th align=center  style="padding:0 15px 0 15px;"> Parameters </th><th> Value </th>'
        for xi in sorted(pvals.keys()):
                html_str += '<tr><td><b>{0}</b></td><td>{1}</td></tr>'.format(
                             xi,
                             pvals[xi])
        html_str += '</table></div>'
        table.value = html_str
        save_button = Button(description='Save Parameter Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.parameter_table.children = [HTML(value='<br>'),
                                         save_button,
                                         table,
                                         make_nominal
                                         ]
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
        table = HTML()
        html_str = '<div><table>\n<caption>Global tolerances determined for Case ' + case.case_number + ' (' + case.signature + ') showing fold-difference to a large qualitative change{0}. </caption>\n'.format(' in log-coordinates' if self.log_coordinates is True else '') 
        html_str += '<tr ><th align=center  rowspan=2 style="padding:0 15px 0 15px;"> Parameters </th><th colspan=2> Tolerance </th></tr>'
        html_str += '<tr><td style="padding:0 15px 0 15px;"><b> Lower bound</b></td><td style="padding:0 15px 0 15px;"><b> Upper bound</b></td></tr>'
        tolerances = case.measure_tolerance(pvals, log_out=self.log_coordinates)
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
        save_button = Button(description='Save Global Tolerance Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.tolerances_table.children = [HTML(value='<br>'),
                                          save_button,
                                          table]
        return
        
        
        
    def identify_parameters(self, b):
        self.pvals = self.case.valid_interior_parameter_set()
        self.update_parameter_table()
        self.update_global_tolerances()
        b.visible = False
                
    def change_logarithmic(self, name, value):
        self.log_coordinates = value
        self.update_equations()
        self.update_global_tolerances()
        
    def change_dynamic_only(self, name, value):
        self.dynamic_only = value
        self.update_equations()
        
    def save_table(self, b):
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        controller.save_widget_data(b)
       
       
        