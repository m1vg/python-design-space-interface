import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg 


class DisplayCase(object):
    
    def __init__(self, controller, case_id, by_signature=False):
        setattr(self, 'controller', controller)
        setattr(self, 'case', controller.ds(case_id, by_signature=by_signature))
        setattr(self, 'info', None)
        setattr(self, 'equations', None)
        setattr(self, 'log_gains', None)
        setattr(self, 'log_coordinates', False)
        
    def create_case_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        self.info = widgets.ContainerWidget()
        self.equations = widgets.ContainerWidget()
        self.log_gains = widgets.ContainerWidget()
        close_button = widgets.ButtonWidget(description='Close Tab')
        close_button.on_click(self.close_widget)
        wi = widgets.PopupWidget(children=[self.info, 
                                               self.equations,
                                               self.log_gains,
                                               close_button])
        wi.set_css('height', '400px')
        self.update_display()
        controller.update_child('Case ' + self.case.case_number, wi)
        
    def close_widget(self, button):
        controller = self.controller 
        controller.update_child('Case ' + self.case.case_number, None)
        
    def update_display(self):
        self.update_info()
        self.update_equations()
        self.update_log_gains()
        
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
        equations = widgets.LatexWidget()
        equations.value = case.equations._repr_latex_()
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
        self.equations.children = [widgets.HTMLWidget(value='<b>Equations:</b><br>'),
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
        html_str = '<br><div><table>\n<caption>Table 2. Logarithmic gains and parameter sensitivities for Case ' + case.case_number + ' (' + case.signature + '). </caption>\n'
        html_str += '<tr ><th rowspan="2" align=center  style="padding:0 15px 0 15px;"> Dependent<br> Variables </th>'
        html_str += '<th colspan="' + str(len(case.independent_variables)) + '" align=center  style="padding:0 15px 0 15px;"> Inependent Variables and Parameters</th></tr><tr align=center>'
        for xi in case.independent_variables:
                html_str += '<td><b>{0}</b></td>'.format(xi)
        html_str += '</tr>\n'
        for xd in case.dependent_variables:
            html_str += '<tr><td align=center  style="padding:0 15px 0 15px;"><b>{0}</b></td>'.format(xd)
            for xi in case.independent_variables:
                html_str += '<td align=center  style="padding:0 15px 0 15px;">{0}</td>'.format(str(case.ssystem.log_gain(xd, xi)))
            html_str += '</tr>\n'
        html_str += '</table></div>'
        table.value = html_str
        self.log_gains.children = [table]
        return
        
    def change_logarithmic(self, name, value):
        self.log_coordinates = value
        self.update_equations()
        
        