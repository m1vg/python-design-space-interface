import dspace
import dspace.plotutils
import dspace.display

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
import mpl_toolkits.mplot3d

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from subprocess import call, Popen, PIPE
from dspace.graphs.designspace_graph import GraphGenerator
from dspace.display.UI.case_widget import DisplayCase
import base64

from collections import OrderedDict

class CaseIntersection(object):
    
    def __init__(self, controller, by_signature=False):
        setattr(self, 'controller', controller)
        setattr(self, 'colocalization_data', VBox())
        setattr(self, 'by_signature', by_signature)
    
    def case_intersection_widget(self):
        controller = self.controller
        cases = Textarea(description='* Cases to intersect:')
        by_signature = Checkbox(description='Cases indicated by signature?', 
                                              value=self.by_signature)
        if 'biological_constraints' not in controller.options:
            bio_constraints = ''
        else:
            bio_constraints = ', '.join(controller.defaults('biological_constraints')) 
        constraints = Textarea(description='Biological constraints:',
                                             value=bio_constraints)
        button = Button(description='Create Case Intersection')
        button.on_click(self.make_case_intersection)
        button.cases = cases
        button.by_signature = by_signature
        button.constraints = constraints
        wi = VBox(children=[cases,
                            by_signature,
                            constraints,
                            button])
        return ('Case Intersections', wi)
        
    def make_case_intersection(self, b):
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
        co = DisplayCaseIntersection(self.controller, cases)
        co.make_display_widget()
        
        
class DisplayCaseIntersection(object):
    
    def __init__(self, controller, cases):
        setattr(self, 'controller', controller)
        setattr(self, 'case_intersection_widget', None)
        setattr(self, 'cases', cases)
        setattr(self, 'name', 'Case Intersection: ' + ', '.join([i.case_number for i in self.cases]))
        setattr(self, 'constraints', [])
        setattr(self, 'active_constraints', False)
        setattr(self, 'ci', dspace.CaseIntersection(cases))
        setattr(self, 'pvals', {})
        
    def make_display_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        self.info = VBox()
        self.constraints_widget = VBox()
        self.plot = HBox()
        self.log_coordinates = True
        self.global_tolerance = VBox()
        close_button = Button(description='Close Tab')
        close_button.on_click(self.close_widget)
        #ss_options = ['log('+ i + ')' for i in controller.ds.dependent_variables]
        #self.y_dropdown = Dropdown(description='y-axis',
        #                           values=ss_options,
        #                           options=ss_options,
        #                           value=self.y_variable)
        #self.make_plot = Button(description='Create Plot')
        #self.make_plot.on_click(self.change_y_axis)
        #self.make_plot.yaxis = self.y_dropdown
        #self.make_plot.visible = True
        check_box = Checkbox(description='Logarithmic coordinates', 
                                           value=self.log_coordinates)
        check_box.on_trait_change(self.update_log, 'value')
        #self.y_dropdown.visible = False
        #if len(self.slice_variables) <= 3:
        #    ## self.make_plot.visible = True
        #    if len(self.slice_variables) == 1:
        #        self.y_dropdown.visible = True
        #    if len(self.slice_variables) == 3:
        #        self.y_dropdown = HTML(value='<font color="red">Warning 3D plots are experimental.</font>')
        wi = VBox(children=[self.info, 
                            self.constraints_widget,
                            check_box,
        #                    self.y_dropdown,
        #                    self.make_plot,
                            self.global_tolerance,
                            close_button])
        if old_ipython is True:
            wi.set_css('height', '400px')
        else:
            wi.height='400px'
            wi.overflow_x = 'auto'
            wi.overflow_y = 'auto'
        self.update_display()
        controller.update_child(self.name, wi)
                
    def update_display(self):
        self.update_info()
        self.update_constraints()
        self.update_global_tolerances()
        
    def update_log(self, name, value):
        controller = self.controller
        self.log_coordinates = value
        self.update_display()
        
        
    def update_info(self):
        
        title = HTML(value='<b> Cases to Co-localize </b>')
        buttons = []
        html_str = '<div><b>Is Valid: {0}</b></div>'.format(self.ci.is_valid())
        valid = HTML(value = html_str)
        pset = self.ci.valid_interior_parameter_set()
        for i in self.cases:
            key = i.case_number
            case_button = Button(description='Case ' + key)
            buttons.append(case_button)
            case_button.pvals = pset[key] if key in pset else None
            case_button.on_click(self.open_case)
        self.info.children = [title] + buttons + [valid]
        
    def update_constraints(self):
        constraints_widget = Textarea(description='Constraints',
                                                    value = ',\n'.join(self.constraints)
                                                    )
        constraints_widget.visible = self.active_constraints
        button = Button(description='Done' if self.active_constraints else 'Modify constraints')
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
        
    
    def update_global_tolerances(self):
        controller = self.controller
        ds = controller.ds
        ci = self.ci
        if ci.is_valid() is False:
            self.global_tolerance.children = []
            return
        pvals = self.ci.valid_interior_parameter_set()
        if pvals is None:
            self.global_tolerance.children = []
            return
        table = HTML()
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
        save_button = Button(description='Save Global Tolerance Table')
        save_button.table_data = html_str
        save_button.on_click(self.save_table)
        self.global_tolerance.children = [HTML(value='<br>'),
                                          save_button,
                                          table]
        return
    
    def save_table(self, b):
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        controller.save_widget_data(b)

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

