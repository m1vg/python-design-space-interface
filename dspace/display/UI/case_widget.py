import dspace
import dspace.plotutils
import dspace.display
from decimal import Decimal
import numpy as np
from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool



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
    from IPython.html.widgets import IntSliderWidget as Slider
    from IPython.html.widgets import FloatSliderWidget as FloatSlider


    VBox = Box
    HBox = Box
    ipy_old = True
else:
    from ipywidgets import *
    from popup import Popup
    ipy_old = False
    
from IPython.display import clear_output, display
from math import log10, floor
from dspace.SWIG.dspace_interface import DSDesignSpaceNumberOfBoundaries

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
        if case_id == '':
            return
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
        if by_signature is True:
            case_id = case_id.strip().replace(' ', '').replace(' ', '')
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
        setattr(self, 'ignore_infinity', True)
        setattr(self, 'tolerance_in_fold_change', True)
        setattr(self, 'dynamic_only', False)
        setattr(self, 'subtitle', subtitle)
        setattr(self, 'custom_parameter_set', False)
        setattr(self, 'check_advanced_settings ', None)
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
        self.eigenvalues_table = VBox()
        self.steady_states_table = VBox()
        self.numerical_boundaries_table = VBox()
        self.calculate_pvals_table = VBox()
        self.bounds_list = []

        self.upper_bound = Slider(description='Upper Boundary Par. Space (Log)',  # 5
                             value=3,
                             min=1,
                             max=20,
                             step=0.25,
                             visible=False)

        self.lower_bound = Slider(description='Lower Boundary Par. Space (Log)',  # 6
                             value=-3,
                             min=-20,
                             max=-1,
                             step=0.25,
                             visible=False)

        self.parameter_method = Dropdown(description='Method',
                                         values=['Tolerances',
                                                 'Vertex Enumeration',
                                                 'Bounding Box',
                                                 'Geometric Mean T. & BB'],
                                         options=['Tolerances',
                                                 'Vertex Enumeration',
                                                 'Bounding Box',
                                                 'Geometric Mean T. & BB'],
                                         value='Tolerances')

        self.enforce_bounds = Checkbox(description='Set Bounds', value=False)
        self.asymmetrical_bounds = Checkbox(description='Asymmetrical Bounds', value=False)
        self.enforce_bounds.on_trait_change(self.update_bounds, 'value')
        self.asymmetrical_bounds.on_trait_change(self.update_bounds, 'value')
        self.bounds_box = Box(children=[self.lower_bound,
                                        self.upper_bound])

        calculate_pvals = Button(description='Determine values for the parameters')
        calculate_pvals.on_click(self.identify_parameters)

        self.calculate_pvals_table.children = [
                                                HTML(value='<br><b> Options for Operating Point Prediction </b>'),
                                               self.parameter_method,
                                               self.enforce_bounds,
                                               self.asymmetrical_bounds,
                                               self.bounds_box,
                                               calculate_pvals,

                                               ]
        self.calculate_pvals_table.visible = False

        if case.is_valid() is True and case.is_cyclical is False:
            # if self.pvals is None:
                self.calculate_pvals_table.visible = True
        close_button = Button(description='Close Tab')
        close_button.on_click(self.close_widget)
        self.check_advanced_settings = Box()
        wi = Popup(children=[self.info,
                             self.equations,
                             self.log_gains,
                             self.bounding_box_table,
                             self.calculate_pvals_table,    #calculate_pvals
                             self.parameter_table,
                             self.tolerances_table,
                             self.check_advanced_settings,
                             self.numerical_boundaries_table,
                             self.eigenvalues_table,
                             self.steady_states_table,
                             close_button]
                   )

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

    def update_bounds(self, value, name):
        if self.enforce_bounds.value is True:
            if self.asymmetrical_bounds.value is False:
                self.lower_bound.visible = True
                self.upper_bound.visible = True
                self.bounds_box.children = [self.lower_bound,
                                           self.upper_bound]
            else:
                ind_variables = self.controller.ds.independent_variables
                asymm_list = []
                children = []
                for variable in ind_variables:
                    asymm_dic = {}
                    asymm_dic['name'] = variable
                    asymm_dic['low_bound'] = FloatSlider(description='Lower Bound (Log)',
                                                    value=-3,
                                                    min=-9,
                                                    max=0,
                                                    step=0.001)

                    asymm_dic['upper_bound'] = FloatSlider(description='Upper Bound (Log)',
                                                      value=3,
                                                      min=0,
                                                      max=9,
                                                      step=0.001,
                                                      )
                    asymm_list.append(asymm_dic)
                    children.append(HTML(value='<b>Bounds for ' + asymm_dic['name'] + '</b>'))
                    children.append(asymm_dic['low_bound'])
                    children.append(asymm_dic['upper_bound'])
                self.bounds_box.children = children
                self.bounds_list = asymm_list
        else:
            self.bounds_box.children = []

        
    def close_widget(self, button):
        controller = self.controller 
        controller.update_child(self.title, None)
        
    def save_parameters(self, b):
        controller = self.controller
        controller.pvals.update(self.pvals)
        self.parameter_table.box_style="success"

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
        html_string += '<b>Case Signature:</b> ' + self.case.signature + '<br>'
        html_string += '<b>Valid:</b> ' + str(self.case.is_valid()) + '<br>'
        html_string += '<b>Is Cyclical:</b> ' + str(self.case.is_cyclical) + '<br>'
        html_string += '<b>Is Blowing:</b> ' + str(self.case.is_unstable) if self.case.is_cyclical is False else '<b>Is Blowing:</b> *'  + '<br><hr>'
        check_box = Checkbox(description='Logarithmic coordinates', 
                                           value=self.log_coordinates)
        check_box.on_trait_change(self.change_logarithmic, 'value')
        html.value = html_string
        self.info.children = [html, check_box]
        
    def update_equations(self):
        # print("processing update_equations")
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
        collapsed_equations = Latex()
        cyclical_description = Latex()
        if self.log_coordinates is True:
            if case.conditions_log is not None:
                conditions_data = case.conditions_log._repr_latex_()
            else:
                conditions_data = None
            solution_data = case.ssystem.solution_log
            boundary_data = case.boundaries_log
        else:
            if case.conditions is not None:
                conditions_data = case.conditions._repr_latex_()
            else:
                conditions_data = None
            solution_data = case.ssystem.solution
            boundary_data = case.boundaries

        if solution_data is not None:
            solution.value = solution_data._repr_latex_()
        else:
            solution.value = 'S-System has no solution.'
        if boundary_data is not None:
            boundaries.value = boundary_data._repr_latex_()
        else:
            boundaries.value = 'S-System has no solution.'
        if conditions_data is not None:
            conditions.value = conditions_data
        else:
            conditions.value = 'S-System has no conditions'

        if case.is_cyclical is True:
            if case.number_of_cycles == 1:
                s = "This cyclical case contains {} cycle. ".format(case.number_of_cycles)
                s += "Its main cyclical variable is {}".format(case.main_cyclical_variables(1).keys()[0]) + " "
                s += "and secondary variable(s): {}.".format(", ".join(case.secondary_cyclical_variables(1).keys()))
            else:
                s = "This cyclical case contains {} cycles. ".format(case.number_of_cycles)
                for i in range(case.number_of_cycles):
                    s += "Cycle {}".format(i+1) + ": \n"
                    s += "Its main cyclical variable is {}".format(case.main_cyclical_variables(i+1).keys()[0]) + ", "
                    s += "its secondary variable(s) {}.".format(", ".join(case.secondary_cyclical_variables(i+1).keys()))
                    s += "\n"
            s += " It has a total of {} subcases, from which {} are valid regular and {} are valid blowing cases.".\
                format(case.number_of_subcases,
                       case.number_of_valid_subcases,
                       case.number_of_valid_blowing_cases)

            cyclical_description.value = s
            collapsed_equations.value = case.augmented_equations._repr_latex_()

        self.equations.children = [check_box,
                                   HTML(value='<b>Equations:</b><br>'),
                                   equations, 
                                   HTML(value='<b>Conditions:</b><br>'),
                                   conditions, 
                                   HTML(value='<b>Solution:</b><br>'),
                                   solution,
                                   HTML(value='<b>Boundaries:</b><br>'),
                                   boundaries] + \
                                   ([
                                    HTML(value='<br>'),
                                    HTML(value='<b>Augmented System:</b><br>'),
                                    HTML(value='<br>'),
                                    collapsed_equations,
                                    cyclical_description,
                                    HTML(value='<br>'),

                                    ] if self.case.is_cyclical is True else [])


    def update_log_gains(self):
        controller = self.controller 
        case = self.case
        if case.ssystem.solution is None: #or case.is_unstable is True:
            self.log_gains.children = []
            return
        table = HTML()

        dependent_variables = case.ssystem.dependent_variables_no_algebraic + case.ssystem.conserved_variables
        html_str = '<div><table>\n<caption>Logarithmic gains and parameter sensitivities for Case ' + case.case_number + ' (' + case.signature + '). </caption>\n'
        html_str += '<tr ><th rowspan="2" align=center  style="padding:0 15px 0 15px;"> Independent Variables <br> and Parameters </th>'
        html_str += '<th colspan="' + str(len(dependent_variables)) + '" align=center  style="padding:0 15px 0 15px;"> Dependent Variables</th></tr><tr align=center>'
        for xd in dependent_variables:
            html_str += '<td><b>{0}</b></td>'.format(xd)
        html_str += '</tr>\n'
        for xi in case.independent_variables:
            html_str += '<tr><td align=center  style="padding:0 15px 0 15px;"><b>{0}</b></td>'.format(xi)
            for xd in dependent_variables:
                html_str += '<td align=center  style="padding:0 15px 0 15px;">{0}</td>'.format(
                    str(case.ssystem.log_gain(xd, xi)))
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

    def string_for_bounding_box(self):

        case = self.case
        tolerances = case.bounding_box(log_out=self.log_coordinates)

        # print("The tolerances are: ", tolerances)

        table = HTML()

        volume = 1
        parameter_count = 0
        max_ratio = None
        min_ratio = []
        neg_vol = False

        lowerBounds = 1e-3
        upperBounds = 1e3

        lower_th = lowerBounds if self.log_coordinates is False else log10(lowerBounds)
        upper_th = upperBounds if self.log_coordinates is False else log10(upperBounds)

        lower_th_inf = 1e-15 if self.log_coordinates is False else -15
        upper_th_inf = 1e15 if self.log_coordinates is False else 15

        html_str = '<div><table>\n<caption>Bounding box determined for Case ' + case.case_number + ' (' + case.signature + '){0}. </caption>\n'.format(' in log-coordinates' if self.log_coordinates is True else '')
        html_str += '<tr ><th align=center  rowspan=2 style="padding:0 15px 0 15px;"> Parameters </th><th colspan=3> Tolerance </th></tr>'
        html_str += '<tr>' \
                    '<td style="padding:0 15px 0 15px;"><b> Lower bound</b></td>' \
                    '<td style="padding:0 15px 0 15px;"><b> Upper bound</b></td>' \
                    '<td style="padding:0 15px 0 15px;"><b> Ratio</b></td>''</tr>'

        for xi in sorted(tolerances.keys()):

            lower, upper = tolerances[xi]
            lower_value = lower if lower > lower_th else lower_th
            upper_value = upper if upper < upper_th else upper_th
            ratio = upper_value / lower_value if self.log_coordinates is False else upper_value - lower_value
            if ratio < 0.0:
                neg_vol = True

            volume = volume * ratio
            parameter_count += 1
            max_ratio = max(max_ratio, ratio)
            min_ratio = min(min_ratio, ratio)

            html_str += '<tr>' \
                        '<td style="padding:0 15px 0 15px;"><b>{0}</b></td>' \
                        '<td style="padding:0 15px 0 15px;">{1}</td>' \
                        '<td style="padding:0 15px 0 15px;">{2}</td>' \
                        '<td style="padding:0 15px 0 15px;">{3}</td>' \
                        '</tr>'.format(
                xi,
                lower if lower > lower_th_inf else '-&infin;',
                upper if upper < upper_th_inf else '&infin;',
                ratio)
        html_str += '</table><caption>'
        html_str += 'Note: Ratios and the bounding box volume are calculated assuming lower and upper bounds for the parameters of ' + '%.1E' % Decimal(
            str(lowerBounds)) + ' and ' + \
                    '%.1E' % Decimal(str(upperBounds)) + ', respectively.'
        try:
            geometric_mean = volume ** (1.0 / parameter_count)
        except:
            geometric_mean = None
        html_str += ' The estimated bounding box volume of the phenotype is ' + '%.2E' % Decimal(str(round(volume,
                                                                                                           3))) + ';' if neg_vol is False else ' The estimated bounding box volume of the phenotype is 0.0;'
        html_str += ' Geometric mean of the ratios: ' + '%.2E' % Decimal(
            str(round(geometric_mean, 3))) + ';' if geometric_mean is not None else ''
        html_str += ' Smallest ratio: ' + '%.2E' % Decimal(
            str(round(min_ratio, 3))) + ';' if min_ratio is not None else ' '
        html_str += ' Largest ratio: ' + '%.2E' % Decimal(str(round(max_ratio, 3))) if max_ratio is not None else ' '
        html_str += '. Ignoring unbounded parameters: False'
        html_str += '. Logarithmic coordinates: ' + str(self.log_coordinates) + '.'
        html_str += '</caption></div>'

        table.value = html_str
        save_button = Button(description='Save Bounding Box Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)

        self.bounding_box_table.children = [HTML(value='<br>'),
                                            save_button,
                                            table,
                                            ]

        return
    
    def update_bounding_box(self):
        controller = self.controller 
        case = self.case
        pvals = self.pvals
        if case.is_valid(p_bounds=pvals) is False:
            self.bounding_box_table.children = []
            return
        if case.is_cyclical is True:
            self.bounding_box_table.children = []
            return
        self.string_for_bounding_box()

        return
        
    def update_parameter_table(self):
        controller = self.controller 
        case = self.case
        pvals = self.pvals
        if pvals is None:
            self.parameter_table.children = []
            return
        if case.is_valid(p_bounds=pvals) is False:
            print("case {} is not valid for the specified parameter range".format(case.case_number))
            self.parameter_table.children = []
            return
        make_nominal = Button(description='Make Nominal Parameter Set')
        make_nominal.on_click(self.save_parameters)
        table = HTML()
        method = self.parameter_method.value + ' method'

        if self.enforce_bounds.value is True:
            if self.asymmetrical_bounds.value is False:
                bounds = ' The bounds in log coordinates are: [' + str(self.lower_bound.value) + ', ' \
                         + str(self.upper_bound.value) + ']'
            else:
                bounds = ' Using asymmetrical bounds.'
        else:
            bounds = 'The bounds in log coordinates are: [-20, 20]'

        html_str = '<div><table>\n<caption>Value for the parameters automatically determined for Case ' + case.case_number + ' (' + case.signature + \
                   ') using the ' + method + '. ' + bounds + '</caption>\n'
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
        pvals = self.pvals if self.custom_parameter_set is False else controller.pvals
        neg_vol = False

        if pvals is None:
            self.tolerances_table.children = []
            return
        if case.is_valid(p_bounds=pvals) is False:
            self.tolerances_table.children = []
            return            
        table = HTML()
        html_str = '<div><table>\n<caption>Global tolerances determined for Case ' + case.case_number + ' (' + case.signature + ') showing fold-difference to a large qualitative change{0}. </caption>\n'.format(' in log-coordinates' if self.log_coordinates is True else '') 
        html_str += '<tr ><th align=center  rowspan=2 style="padding:0 15px 0 15px;"> Parameters </th><th colspan=3> Tolerance </th></tr>'
        html_str += '<tr>' \
                    '<td style="padding:0 15px 0 15px;"><b> Lower bound</b></td>' \
                    '<td style="padding:0 15px 0 15px;"><b> Upper bound</b></td>' \
                    '<td style="padding:0 15px 0 15px;"><b> Ratio</b></td>''</tr>'
        tolerances = case.measure_tolerance(pvals, log_out=self.log_coordinates)
        volume = 1
        parameter_count = 0
        max_ratio = None
        min_ratio = []

        lowerBounds = 1e-3
        upperBounds = 1e3

        lower_th = lowerBounds if self.log_coordinates is False else log10(lowerBounds)  #1e-15
        upper_th = upperBounds if self.log_coordinates is False else log10(upperBounds)   #1e15

        lower_th_inf = 1e-15 if self.log_coordinates is False else -15
        upper_th_inf = 1e15 if self.log_coordinates is False else 15

        for xi in sorted(pvals.keys()):

            lower, upper = tolerances[xi]
            lower = lower * pvals[xi] if self.log_coordinates is False else lower + log10(pvals[xi])
            upper = upper * pvals[xi] if self.log_coordinates is False else upper + log10(pvals[xi])

            lower_value = lower if lower > lower_th else lower_th
            upper_value = upper if upper < upper_th else upper_th
            ratio = upper_value/lower_value if self.log_coordinates is False else upper_value - lower_value

            if ratio < 0.0:
                neg_vol = True

            if self.ignore_infinity is True:
                if lower > lower_th and upper < upper_th:
                    volume = volume*ratio
                    parameter_count += 1
                    max_ratio = max(max_ratio, ratio)
                    min_ratio = min(min_ratio, ratio)
            else:
                volume = volume * ratio
                parameter_count += 1
                max_ratio = max(max_ratio, ratio)
                min_ratio = min(min_ratio, ratio)

            lower, upper = tolerances[xi]
            if self.tolerance_in_fold_change is False:
                lower = lower * pvals[xi] if self.log_coordinates is False else lower + log10(pvals[xi])
                upper = upper * pvals[xi] if self.log_coordinates is False else upper + log10(pvals[xi])
            html_str += '<tr>' \
                        '<td style="padding:0 15px 0 15px;"><b>{0}</b></td>' \
                        '<td style="padding:0 15px 0 15px;">{1}</td>' \
                        '<td style="padding:0 15px 0 15px;">{2}</td>' \
                        '<td style="padding:0 15px 0 15px;">{3}</td>' \
                        '</tr>'.format(
                         xi,
                         lower if lower > lower_th_inf else '-&infin;',
                         upper if upper < upper_th_inf else '&infin;',
                         ratio)

        html_str += '</table><caption>'
        html_str += 'Note: Global tolerance calculated based on the following values for the parameters: ' 
        html_str += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals.keys())]) + '.' \
                    ' Ratios and the tolerance volume are calculated assuming lower and upper bounds for the parameters of ' + '%.1E' % Decimal('0.001') + ' and '+ \
                    '%.1E' % Decimal('1000') + ', respectively.'
        if parameter_count == 0:
            html_str += ' The volume of the phenotype and the geometric mean of the ratios are not defined since all parameters are unbounded.'
        else:
            try:
                geometric_mean = volume**(1.0/parameter_count)
            except:
                geometric_mean = None
            html_str += ' The estimated tolerance volume of the phenotype is ' + '%.2E' % Decimal(str(round(volume, 3))) + ';' if neg_vol is False else ' The estimated volume of the phenotype is 0.0;'
            html_str += ' Geometric mean of the ratios: ' + '%.2E' % Decimal(str(round(geometric_mean, 3))) + ';' if geometric_mean is not None else ''
            html_str += ' Smallest ratio: ' + '%.2E' % Decimal(str(round(min_ratio, 3))) + ';' if min_ratio is not None else ' '
            html_str += ' Largest ratio: ' + '%.2E' % Decimal(str(round(max_ratio, 3)))  if max_ratio is not None else ' '
            html_str += '. Ignoring unbounded parameters: ' + str(self.ignore_infinity) + '.' #if self.log_coordinates is False else ' '

        html_str += '</caption></div>'
        table.value = html_str
        save_button = Button(description='Save Global Tolerance Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.tolerances_table.children = [HTML(value='<br>'),
                                          save_button,
                                          table,
                                          ]
        return

    def update_eigen_values(self):

        controller = self.controller
        case = self.case
        pvals = self.pvals if self.custom_parameter_set is False else controller.pvals

        try:
            ssys = case.ssystem.remove_algebraic_constraints()
            eigen_values = ssys.eigenvalues(pvals)
            sort_index = np.argsort(eigen_values.real)
            num_pos = len([x for x in eigen_values.real if x > 0])

            # print("The real part of the eigenvalues is: ", eigen_values.real)
            # print("The imaginary part of the eigenvalues is: ", eigen_values.imag)

            if pvals is None:
                self.eigenvalues_table.children = []
                return
            if case.is_valid(p_bounds=pvals) is False:
                self.eigenvalues_table.children = []
                return

            table = HTML()
            html_str = '<div><caption>Eigenvalues determined for Case ' + case.case_number + ' (' + case.signature + ') {0}. </caption><table>'.format(' in log-coordinates' if self.log_coordinates is True else '')
            html_str += '<tr ><th align=center  rowspan=2 style="padding:0 15px 0 15px;"> Rank </th><th colspan=2> Eigenvalues </th></tr>'
            html_str += '<tr>' \
                        '<td style="padding:0 15px 0 15px;"><b> Real</b></td>' \
                        '<td style="padding:0 15px 0 15px;"><b> Imaginary</b></td>''</tr>'
            count = 0
            for x in sort_index:
                html_str += '<tr>' \
                            '<td style="padding:0 15px 0 15px;"><b>{0}</b></td>' \
                            '<td style="padding:0 15px 0 15px;">{1}</td>' \
                            '<td style="padding:0 15px 0 15px;">{2}</td>' \
                            '</tr>'.format(
                             count+1,
                             eigen_values.real[x],
                             eigen_values.imag[x])
                count += 1
            html_str += '</table><caption>'
            html_str += 'Note: Eigenvalues calculated based on the following values for the parameters: '
            html_str += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals.keys())]) + '.'
            html_str += ' The number of eigenvalues with positive real part is {0}. '.format(num_pos)
            html_str += '</caption></div>'
            table.value = html_str
            save_button = Button(description='Save Eigenvalues Table')
            save_button.table_data = html_str
            save_button.on_click(self.save_table)
            self.eigenvalues_table.children = [HTML(value='<br>'),
                                              save_button,
                                              table,
                                              ]
        except:
            self.eigenvalues_table.children = []
            return

    def update_steady_states(self):

        controller = self.controller
        case = self.case
        pvals = self.pvals if self.custom_parameter_set is False else controller.pvals

        try:
            if pvals is None:
                self.steady_states_table.children = []
                return
            if case.is_valid(p_bounds=pvals) is False:
                self.steady_states_table.children = []
                return
            steady_states = case.steady_state(pvals)
            table = HTML()
            html_str = '<div><caption>Steady States determined for Case ' + case.case_number + ' (' + case.signature + ') {0}. </caption><table>'.format(
                ' in log-coordinates' if self.log_coordinates is True else '')
            html_str += '<tr ><th align=center  style="padding:0 15px 0 15px;"> Variable </th><th> SS Value </th>'
            for xi in sorted(steady_states.keys()):
                html_str += '<tr><td><b>{0}</b></td><td>{1}</td></tr>'.format(
                    xi,
                    steady_states[xi] if self.log_coordinates is False else np.log10(steady_states[xi]))
            html_str += '</table></div>'
            html_str += '</table><caption>'
            html_str += 'Note: Steady States calculated based on the following values for the parameters: '
            html_str += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals.keys())]) + '.'
            html_str += '</caption></div>'
            table.value = html_str
            save_button = Button(description='Save Steady States Table')
            save_button.table_data = html_str
            save_button.on_click(self.save_table)
            self.steady_states_table.children = [HTML(value='<br>'),
                                                 save_button,
                                                 table,
                                                 ]
        except:
            self.steady_states_table.children = []
            return

    def update_double_value_boundaries_at_point(self):

        controller = self.controller
        case = self.case
        pvals = self.pvals if self.custom_parameter_set is False else controller.pvals
        log_out = self.log_coordinates
        if pvals is None:
            return
        if case.is_valid(p_bounds=pvals) is False:
            return

        boundaries_at_point = case.double_value_boundaries_at_point(pvals, log_out=log_out)
        table = HTML()
        html_str = '<div><caption>Numerical Boundaries for Case ' + case.case_number + ' (' + case.signature + ') {0}. </caption><table>'.format(
            ' in log-coordinates' if self.log_coordinates is True else '')
        html_str += '<tr ><th align=center  style="padding:0 15px 0 15px;"> Boundary </th><th> Value </th>'
        on_edge = []
        for index, value in enumerate(boundaries_at_point):
            html_str += '<tr><td><b>{0}</b></td><td>{1}</td></tr>'.format(index + 1,
                                                                          value)
            if log_out is False:
                if round(value) == 1.0:
                    on_edge.append(str(index + 1))
            else:
                if value < 1e-14:
                    on_edge.append(str(index + 1))
        html_str += '</table></div>'
        html_str += '</table><caption>'
        html_str += 'Note: Numerical values for boundaries calculated based on the following values for the parameters: '
        html_str += '; '.join([i + ' = ' + str(pvals[i]) for i in sorted(pvals.keys())]) + '.'
        if len(on_edge) != 0:
            html_str += ' Following boundaries are on edge: ' + ', '.join(on_edge)
        html_str += '</caption></div>'
        table.value = html_str
        save_button = Button(description='Save Numerical Boundaries Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.numerical_boundaries_table.children = [HTML(value='<br>'),
                                                    save_button,
                                                    table,
                                                    ]
        
    def ignore_unbounded_parameters(self,name,value):
        self.ignore_infinity = value
        try:
            self.update_global_tolerances()
        except:
            pass

    def update_tolerance_table(self, name, value):
        self.tolerance_in_fold_change = value
        try:
            self.update_global_tolerances()
        except:
            pass

    def use_custom_parameter_set(self, name, value):
        self.ignore_infinity = self.check_advanced_settings.children[1].value
        self.custom_parameter_set = value
        try:
            self.update_global_tolerances()
        except:
            pass
        try:
            self.update_eigen_values()
        except:
            pass
        try:
            self.update_steady_states()
        except:
            pass
        try:
            self.update_double_value_boundaries_at_point()
        except:
            pass

    def process_lrs_volume(self, b):

        log_coordinates = self.log_coordinates

        if b.description == 'Expand lrs Options':
            b.description = 'Calculate lrs Volume'
            self.check_advanced_settings.upper_bound.visible = True
            self.check_advanced_settings.lower_bound.visible = True
            self.check_advanced_settings.nr_vertices_checkbox.visible = True
        elif b.description == 'Calculate lrs Volume':
            self.check_advanced_settings.html_vol_table.visible = True
            limit_vertices = self.check_advanced_settings.limit_vertices_checkbox.value
            lb = self.controller.pvals.copy()
            ub = self.controller.pvals.copy()
            maxVertices = int(self.check_advanced_settings.max_nr_vertices.value)
            for key in lb.keys():
                lb[key] = 10**self.check_advanced_settings.lower_bound.value
                ub[key] = 10**self.check_advanced_settings.upper_bound.value

            if limit_vertices is True:
                limiting = 'Max. {} vertices: '.format(self.check_advanced_settings.max_nr_vertices.value)
            else:
                limiting = ''
            volume, nr_vertices, vertices_matrix, operating_point = self.case.volume_lrs(lb, ub, maxVertices,
                                                                                         limit_vertices,
                                                                                         return_vertices_matrix=True)
            print("The volume is: ", volume)
            volume = volume if nr_vertices != 0 else 0
            volume = '%.2E' % Decimal(str(round_sig(volume, 3))) if volume != 0 else '0.0'


            html_string = '<div><table>\n<caption>Vertices for Case ' + self.case.case_number + ' (' + self.case.signature + '). </caption>\n'
            html_string += '<tr ><th rowspan="2" align=center  style="padding:0 15px 0 15px;"> Independent Variables <br> and Parameters </th>'

            html_string += '<th colspan="' + str(nr_vertices) + '" align=center  style="padding:0 15px 0 15px;"> Vertices '

            if log_coordinates is True:
                html_string +='(Log) </th></tr><tr align=center>'
            else:
                html_string += '</th></tr><tr align=center>'

            for vertex in range(int(nr_vertices)):
                html_string += '<td><b>{0}</b></td>'.format(vertex+1)
            html_string += '<td><b>{0}</b></td>'.format("Average<br>(Calculated in log coordinates)")
            html_string += '</tr>\n'
            for var_nr, xi in enumerate(self.case.independent_variables):
                html_string += '<tr><td align=center  style="padding:0 15px 0 15px;"><b>{0}</b></td>'.format(xi)
                for vertex in range(int(nr_vertices)):
                    if log_coordinates is True:
                        value = vertices_matrix[vertex][var_nr]
                    else:
                        value = 10**vertices_matrix[vertex][var_nr]
                    html_string += '<td align=center  style="padding:0 15px 0 15px;">{0}</td>'.format(str(value))
                if log_coordinates is True:
                    value = operating_point[xi]
                else:
                    value = 10**operating_point[xi]
                html_string += '<td align=center  style="padding:0 15px 0 15px;">{0}</td>'.format(str(value))
                html_string += '</tr>\n'
            html_string += '</table></div><br>'
            html_string += limiting + 'The phenotype had {} vertices and a lrs log volume of {}. <br>'.format(nr_vertices,
                                                                                                                 volume)

            save_button = Button(description='Save Vertices Table')
            save_button.table_data = html_string
            save_button.on_click(self.save_table)
            self.check_advanced_settings.html_vol_table.children = [
                                                                    save_button,
                                                                    HTML(value=html_string)
                                                                    ]

    def update_settings_box(self):
        if self.check_advanced_settings.nr_vertices_checkbox.value is True:
            self.check_advanced_settings.max_nr_vertices.visible = True
        else:
            self.check_advanced_settings.max_nr_vertices.visible = False

    def update_check_advanced_settings_box(self):
        controller = self.controller
        nr_ineq = DSDesignSpaceNumberOfBoundaries(controller.ds._swigwrapper) + \
                  len(controller.constraints) + len(controller.pvals)*2
        max_vertices = nr_ineq ** (len(controller.pvals) / 2)

        check_box0 = Checkbox(description='Report Tolerance in Fold-change', value=True)
        check_box0.on_trait_change(self.update_tolerance_table, 'value')

        check_box1 = Checkbox(description='Ignore unbounded parameters for volume determination:',
                              value=True)
        check_box1.on_trait_change(self.ignore_unbounded_parameters, 'value')


        check_box2 = Checkbox(description='Use customized parameter set:',
                              value=False)
        check_box2.on_trait_change(self.use_custom_parameter_set, 'value')


        limit_vertices_checkbox = Checkbox(description='Limit Number of Vertices',
                                        value=False,
                                        visible=False)

        max_nr_vertices = Slider(description='Max. Number of Vertices',  # 4
                                 value=max_vertices,
                                 min=1,
                                 max=round_sig(max_vertices, 1)*2,
                                 step=round_sig(max_vertices, 1)/100,
                                 visible=False)

        upper_bound = Slider(description='Upper Boundary Par. Space (Log)',  # 5
                             value=4,
                             min=1,
                             max=20,
                             step=0.25,
                             visible=False)

        lower_bound = Slider(description='Lower Boundary Par. Space (Log)',  # 6
                             value=-4,
                             min=-20,
                             max=-1,
                             step=0.25,
                             visible=False)

        html_vol_table = VBox()

        get_vertices_button = Button(description='Expand lrs Options', visible=True)
        get_vertices_button.on_click(self.process_lrs_volume)
        limit_vertices_checkbox.on_trait_change(self.update_settings_box, 'value')

        self.check_advanced_settings.check_box1 = check_box1
        self.check_advanced_settings.limit_vertices_checkbox = limit_vertices_checkbox
        self.check_advanced_settings.check_box2 = check_box2
        self.check_advanced_settings.upper_bound = upper_bound
        self.check_advanced_settings.lower_bound = lower_bound
        self.check_advanced_settings.max_nr_vertices = max_nr_vertices
        self.check_advanced_settings.nr_vertices_checkbox = limit_vertices_checkbox
        self.check_advanced_settings.get_vertices_button = get_vertices_button
        self.check_advanced_settings.html_vol_table = html_vol_table

        if self.case.is_valid(p_bounds=self.controller.pvals) is True:
            self.check_advanced_settings.children = [check_box0,
                                                     check_box1,
                                                     check_box2,
                                                     get_vertices_button,
                                                     lower_bound,
                                                     upper_bound,
                                                     limit_vertices_checkbox,
                                                     max_nr_vertices,
                                                     html_vol_table
                                                     ]
        else:
            self.check_advanced_settings.children = [check_box0,
                                                     check_box1,
                                                     get_vertices_button,
                                                     lower_bound,
                                                     upper_bound,
                                                     limit_vertices_checkbox,
                                                     max_nr_vertices,
                                                     html_vol_table
                                                     ]

    def identify_parameters(self, b):

        controller = self.controller
        if self.enforce_bounds.value is True:
            var = VariablePool(names=controller.ds.independent_variables)
            p_bounds = {}
            if self.asymmetrical_bounds.value is False:
                for i in var:
                    p_bounds[i] = [10**self.lower_bound.value, 10**self.upper_bound.value]
            else:
                for var_dic in self.bounds_list:
                    p_bounds[var_dic['name']] = [10**(var_dic['low_bound'].value), 10**(var_dic['upper_bound'].value)]
        else:
            p_bounds = None

        if self.parameter_method.value == 'Tolerances':
            try:
                self.pvals = self.case.valid_interior_parameter_set(p_bounds=p_bounds)
            except:
                print("Case {} is not valid within the specified parameter range. Try [-20, 20]".format(self.case.case_number))
        elif self.parameter_method.value == 'Vertex Enumeration':
            self.pvals = self.case.valid_interior_parameter_set_vertex_enumeration(p_bounds=p_bounds)
        elif self.parameter_method.value == 'Bounding Box':
            self.pvals = self.case.valid_interior_parameter_bounding_box(p_bounds=p_bounds)
        elif self.parameter_method.value == 'Geometric Mean T. & BB':
            self.pvals = self.case.valid_interior_parameter_geometric_mean_t_bb(p_bounds=p_bounds)

        if self.pvals is None or self.case.is_valid(p_bounds=p_bounds) is False:
            print("case {} is not valid for the specified parameter range".format(self.case.case_number))
            return

        self.update_parameter_table()
        self.update_global_tolerances()
        self.update_eigen_values()
        self.update_steady_states()
        self.update_double_value_boundaries_at_point()
        # b.visible = False
        self.update_check_advanced_settings_box()

    def change_logarithmic(self, name, value):
        self.log_coordinates = value
        self.update_equations()
        self.update_bounding_box()
        try:
            b = VBox()
            b.description = 'Calculate lrs Volume'
            self.process_lrs_volume(b)
        except:
            pass
        try:
            self.update_global_tolerances()
        except:
            pass
        self.update_steady_states()
        self.update_double_value_boundaries_at_point()

    def change_dynamic_only(self, name, value):
        self.dynamic_only = value
        self.update_equations()
        
    def save_table(self, b):
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        controller.save_widget_data(b)


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)
