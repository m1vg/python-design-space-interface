import dspace
from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
import dspace.plotutils
import dspace.display
# from InteractiveInput import load_widget
import cPickle as pickle
from dspace.expressions import Expression
import numpy as np
from scipy.integrate import odeint
from distutils.version import LooseVersion, StrictVersion
import IPython
from math import *
import pandas as pd
import difflib
import random
np.seterr(all='raise')
from dspace.display.UI.detect_peaks import detect_peaks

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
    from IPython.html.widgets import FloatSliderWidget as FloatSlider
    VBox = Box
    HBox = Box
    old_ipython = True
else:
    from ipywidgets import *
    from popup import Popup
    old_ipython = False
    
from IPython.display import clear_output, display
import matplotlib.pyplot as plt

import cStringIO
import base64
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from subprocess import call, Popen, PIPE
from dspace.graphs.designspace_graph import GraphGenerator
import itertools

# This is to import assimulo without warning regarding missing solver that won't be used anyway.
from IPython.utils import io
with io.capture_output() as capture:
    from assimulo.solvers import IDA
    from assimulo.solvers import Radau5DAE
    from assimulo.problem import Implicit_Problem



class FullSystem(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'plot_data', VBox())
        setattr(self, 'title', None)
        setattr(self, 'caption', None)
        if 'full_system_options' in controller.options.keys():
                options_dic = controller.options['full_system_options']
                for key in options_dic.keys():
                    setattr(self, key, options_dic[key])
                keys = ['response_times', 'dynamical_phenotypes', 'phase_shifts']
                for key in keys:
                    if key not in options_dic.keys():
                        setattr(self, key, True)
                if 'initial_conditions' not in options_dic.keys():
                    setattr(self, 'initial_conditions', None)
        else:
                keys = ['response_times', 'dynamical_phenotypes', 'phase_shifts']
                for key in keys:
                    setattr(self, key, True)
                setattr(self, 'initial_conditions', None)

        self.coef_dic = {}
        self.conservations_dictionary = None

        # 1. Get the auxiliary, independent and dependent variables
        auxiliary = controller.ds.auxiliary_variables
        dependent = controller.ds.dependent_variables
        dependent_no_artificial = controller.ds.dependent_variables

        # 2. Remove auxiliary variables from dependent variables
        for aux in auxiliary:
            if aux in dependent:
                dependent.remove(aux)

        # 3. Generate conserved, artificial variables Xc1, ... Xcn
        if self.controller.ds.number_conservations != 0:
            artificial_var = ['Xc' + str(n+1) for n in range(int(self.controller.ds.number_conservations))]

            # Remove artificial variables from dependent variables
            for aux in artificial_var:
                if aux in dependent_no_artificial:
                    dependent_no_artificial.remove(aux)

        setattr(self, 'dependent_no_aux', dependent)
        self.generate_ode_equations()
        dependent_no_aux_no_conserved = [variable for variable in dependent if variable not in self.coef_dic.keys()]
        setattr(self, 'dependent_no_aux_no_conserved', dependent_no_aux_no_conserved)
        self.dependent_no_artificial = dependent_no_artificial

    def show_dae_system(self, b):

        # This function generates a pop-up element

        close_button = Button(description='Close')
        alge_string = [str(i) for i in self.conservations]
        ode_string = [str(i) for i in self.ode_equations]
        eqs_string = ode_string + alge_string
        cons_latex = Latex(value=dspace.Equations(eqs_string, latex_symbols=self.controller.symbols)._repr_latex_())
        content = Box(children=[cons_latex,
                                close_button])
        p = Popup(children=[content], description='System of Differential Algebraic Equations')
        display(p)
        close_button.window = p
        close_button.on_click(lambda x: x.window.close())
        
    def fullsystem_widget(self):
        controller = self.controller
        plot_type = Dropdown(description='Plot Type',
                             values=['Time Courses', 'Titrations', 'Trajectories', 'Event Integration'],
                             options=['Time Courses', 'Titrations', 'Trajectories', 'Event Integration'],
                             value='Time Courses')

        x_axis = Dropdown(description='X-Axis',
                          values=['Time'],
                          options=['Time'],
                          value='Time')

        y_axis = Dropdown(description='Y-Axis',
                          values=self.dependent_no_artificial, # was self.dependent_no_aux
                          options=self.dependent_no_artificial,
                          value=self.dependent_no_artificial[0])

        y_axis_secondary = Dropdown(description='Secondary Y-Axis',
                                    values=self.dependent_no_artificial + ['None'],
                                    options=self.dependent_no_artificial + ['None'],
                                    value='None')

        self.solver = Dropdown(description='DAE Solver',
                               values=['IDA', 'RADAU5'],
                               options=['IDA', 'RADAU5'],
                               value='IDA')

        self.automatic_scale = True
        automatic_scale = Button(description='Set Scales')
        self.automatic_scale_button = automatic_scale
        self.automatic_scale_button.on_click(self.update_scales)
        self.refresh_initial_conditions = Button(description='Refresh Initial Conditions')
        self.refresh_initial_conditions.on_click(self.update_initial_conditions)

        self.coordinates = Dropdown(description='Coordinates',
                                    values=['Cartesian', 'Logarithmic'],
                                    options=['Cartesian', 'Logarithmic'],
                                    value='Logarithmic'
                                    )

        x_min = FloatText(description='X-Min', value=0.01, visible=False)
        x_max = FloatText(description='X-Max', value=100, visible=False)

        y_min = FloatText(description='log Y-Min', value=0.01, visible=False)
        y_max = FloatText(description='log Y-Max', value=100, visible=False)

        y_min_secondary = FloatText(description='Secondary Y-Min', value=0.01, visible=False)
        y_max_secondary = FloatText(description='Secondary Y-Max', value=100, visible=False)

        x_min_time = FloatText(description='Initial Time Point', value=0, visible=True)
        x_max_time = FloatText(description='Final Time Point', value=100, visible=True)
        self.step_time = FloatText(value=1000, description='Number of Points for Time')

        title_s = 'Analysis of the ' + controller.name + ' by system design space showing a time course plot'
        title = Text(description='Title', value=title_s)
        caption = Textarea(description='Caption', value='Concentration of ' + y_axis.value + ' is plotted against time.')

        add_plot = Button(description='Add Plot')
        add_plot.on_click(self.call_ode_case)

        show_dae = Button(description='Show DAE System', visible=False)
        show_dae.on_click(self.show_dae_system)
        if self.controller.conserved is True:
            show_dae.visible = True

        full_widget = Box()
        advanced_box = Box()
        full_widget.children = [plot_type, x_axis, y_axis, y_axis_secondary, x_min, x_max,
                                y_min, y_max, y_min_secondary, y_max_secondary, x_min_time,
                                x_max_time, self.step_time,
                                advanced_box, title, caption, add_plot, show_dae]

        # Attach all elements to the class.
        self.type = plot_type
        self.xlabel = x_axis
        self.ylabel = y_axis
        self.ylabel_secondary = y_axis_secondary
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.y_min_secondary = y_min_secondary
        self.y_max_secondary = y_max_secondary
        self.x_min_time = x_min_time
        self.x_max_time = x_max_time


        self.title = title
        self.caption = caption
        self.advanced_box = advanced_box

        # we create an advanced box for each of the three cases and attach it to the main advanced box. .visible=False
        # is used to control which box is shown.

        # Time courses
        self.advanced_time_courses = Box(visible=True)
        initial_string = HTML(value='<b>Initial Concentrations</b>')
        self.initial_concentration = dict((i, FloatText(description=i, value=0.001)) for i in self.dependent_no_aux)
        self.time_courses_copy_final_box = Box()

        response_time = Text(value='', description='Calculate Response Time for Variable or Function')
        response_time_title = HTML(value='<b>Response Times</b>')
        response_time_method = Dropdown(description='Method for Target y-value',
                                        values=['y_final - y_initial',
                                                'max(y) - min(y)',
                                                'y_final*(1-Thr.) ; y_final*(1+Thr.)',
                                                'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ',
                                                'custom'],
                                        options=['y_final - y_initial',
                                                 'max(y) - min(y)',
                                                 'y_final*(1-Thr.) ; y_final*(1+Thr.)',
                                                 'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ',
                                                 'custom'],
                                        value='y_final - y_initial'
                                        )

        ignore_if_less = FloatSlider(description='Threshold value (%) :',
                                     value=5, min=0, max=100, step=0.1, visible=False
                                     )
        update_if_less = lambda name, value: update_if_less_function(name, value, widget=ignore_if_less)
        response_time_method.on_trait_change(update_if_less, 'value')

        # fun_update = lambda name, value: self.update_visibility(name, value, widget=self.custom_points_titration)
        # self.checkbox_custom_points_titration.on_trait_change(fun_update, 'value')

        response_time_box = Box(children=[response_time_title,
                                          response_time,
                                          response_time_method,
                                          ignore_if_less
                                          ],
                                visible=self.response_times)

        add_ss_function_button = Button(description='Add Steady State Function')
        ss_functions = Box()
        ss_functions.members = []

        self.advanced_time_courses.response_time = response_time
        self.advanced_time_courses.response_time_method = response_time_method
        self.advanced_time_courses.ignore_if_less = ignore_if_less
        self.advanced_time_courses.add_ss_function_button = add_ss_function_button
        self.advanced_time_courses.ss_functions = ss_functions

        add_ss_function_button.on_click(self.populate_ss_functions)

        self.dynamical_phenotypes_box = Box(visible=self.dynamical_phenotypes)
        dynamical_phenotypes_string = HTML(value='<b>Dynamical Phenotypes</b>')
        self.include_dynamical_in_excel = Checkbox(value=False, description='Include Dynamical Phenotypes in .xlsx File')
        self.include_dynamical_in_plot = Checkbox(value=False, description='Show Dynamical Phenotypes in Plot')

        self.dynamical_phenotypes_box.children = [dynamical_phenotypes_string,
                                                  self.include_dynamical_in_excel,
                                                  self.include_dynamical_in_plot
                                                  ]

        self.phase_shift_box = Box(visible=self.phase_shifts)
        phase_shift_string = HTML(value='<b>Phase Shifts</b>')
        self.reference_time_course = Dropdown(values=['None'] + self.dependent_no_artificial,
                                              options=['None'] + self.dependent_no_artificial,
                                              value='None',
                                              description='Reference Time Course')

        self.target_phase_shift_course = Dropdown(values=['None'] + self.dependent_no_artificial,
                                                  options=['None'] + self.dependent_no_artificial,
                                                  value='None',
                                                  description='Calculate Phase Shift for Time Course')

        self.phase_shift_box.children = [phase_shift_string,
                                         self.reference_time_course,
                                         self.target_phase_shift_course]

        if controller.ds.number_conservations != 0:
            self.advanced_time_courses.children = [initial_string] + \
                                                  [value for key, value in self.initial_concentration.iteritems()] + \
                                                  [self.time_courses_copy_final_box] + \
                                                  [self.refresh_initial_conditions] +\
                                                  [self.automatic_scale_button] + \
                                                  [self.coordinates] +\
                                                  [ss_functions] + \
                                                  [add_ss_function_button] + \
                                                  [response_time_box] + \
                                                  [self.solver] +\
                                                  [self.dynamical_phenotypes_box] +\
                                                  [self.phase_shift_box]

        else:
            self.advanced_time_courses.children = [initial_string] + \
                                                  [value for key, value in self.initial_concentration.iteritems()] + \
                                                  [self.time_courses_copy_final_box] + \
                                                  [self.refresh_initial_conditions] + \
                                                  [self.automatic_scale_button] + \
                                                  [self.coordinates] + \
                                                  [ss_functions] + \
                                                  [add_ss_function_button] + \
                                                  [response_time_box] +\
                                                  [self.dynamical_phenotypes_box] +\
                                                  [self.phase_shift_box]

        # Titrations
        self.step= FloatText(value=100, description='Number of Points for ' + str(self.xlabel.value))
        self.backwards = Checkbox(value=False, description='Show Backwards Integration')
        self.s_system_steady_states_box = Box()
        self.s_system_steady_states_box.show_ssystems = Checkbox(value=False, description='Show S-System Approximation')
        self.s_system_steady_states_box.only_cases = Textarea(description='Only Cases', visible=False)
        fun_update2 = lambda name, value: self.update_visibility2(name, value, widget=self.s_system_steady_states_box.only_cases)
        self.s_system_steady_states_box.show_ssystems.on_trait_change(fun_update2, 'value')
        self.s_system_steady_states_box.children = [self.s_system_steady_states_box.show_ssystems,
                                                    self.s_system_steady_states_box.only_cases,
                                                    ]
        self.titration_function = Text(value='', description='Custom Function', visible=False)

        self.checkbox_custom_points_titration = Checkbox(value=False, description='Custom Values for X-axis')
        self.custom_points_titration = Textarea(description='Values', value='1, 1.00006244, 1.000124883, 1.00018733, 1.000249782, 1.000312237, 1.000374696, 1.000437159, 1.000499625, 1.000562096, 1.000624571, 1.000687049, 1.000749532, 1.000812018, 1.000874508, 1.000937002, 1.0009995, 1.001062002, 1.001124508, 1.001187018, 1.001249532, 1.001312049, 1.001374571, 1.001437096, 1.001499625, 1.001562158, 1.001624696, 1.001687237, 1.001749781, 1.00181233, 1.001874883, 1.001937439, 1.002, 1.002849895, 1.003700511, 1.004551849, 1.005403908, 1.006256691, 1.007110197, 1.007964426, 1.00881938, 1.00967506, 1.010531465, 1.011388596, 1.012246455, 1.013105041, 1.013964356, 1.014824399, 1.015685172, 1.016546675, 1.017408908, 1.018271873, 1.01913557, 1.02, 1.036712367, 1.05369856, 1.070963067, 1.088510447, 1.106345335, 1.124472442, 1.142896555, 1.161622542, 1.180655348, 1.2',
                                                visible=False)
        fun_update = lambda name, value: self.update_visibility(name, value, widget=self.custom_points_titration)
        self.checkbox_custom_points_titration.on_trait_change(fun_update, 'value')

        self.advanced_titrations = Box(visible=False)
        self.initial_concentration2 = dict((i, FloatText(description=i, value=0.001)) for i in self.dependent_no_aux)

        self.update_initial_conditions(None)

        export_excel_data = Button(visible=False, description='Export Results to .xlsx File')
        self.advanced_titrations.export_excel_data = export_excel_data

        self.response_times_box = Box()
        calculate_response_time_titration = Checkbox(description='Calculate Response Times', value=False)
        calculate_response_time_titration.on_trait_change(self.update_response_time_box, 'value')
        response_time_method_titration = Dropdown(description='Method for Target y-value',
                                                  values=['y_final - y_initial',
                                                          'max(y) - min(y)',
                                                          'y_final*(1-Thr.) ; y_final*(1+Thr.)',
                                                          'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ',
                                                          'custom'],
                                                  options=['y_final - y_initial',
                                                           'max(y) - min(y)',
                                                           'y_final*(1-Thr.) ; y_final*(1+Thr.)',
                                                           'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ',
                                                           'custom'],
                                                  value='y_final - y_initial',
                                                  visible=False
                                                )

        response_time_ignore_if_less_titration = FloatSlider(description='Threshold value (%) :',
                                                             value=5,
                                                             min=0,
                                                             max=100,
                                                             step=0.1,
                                                             visible=False
                                                            )

        update_if_less = lambda name, value: update_if_less_function(name, value,
                                                                     widget=response_time_ignore_if_less_titration)
        response_time_method_titration.on_trait_change(update_if_less, 'value')

        reference_state = Dropdown(description='Reference Point for Response Time',
                                   values=['Initial Point', 'Previous Point'],
                                   options=['Initial Point', 'Previous Point'],
                                   value='Initial Point',
                                   visible=False)

        self.response_times_box.children = [calculate_response_time_titration,
                                            response_time_method_titration,
                                            response_time_ignore_if_less_titration,
                                            reference_state]

        self.phase_shift_box2 = Box()
        calculate_phase_shift_titration = Checkbox(description='Calculate Phase Shift', value=False)
        calculate_phase_shift_titration.on_trait_change(self.update_phase_shift_box, 'value')
        reference_time_course = Dropdown(values=['None'] + self.dependent_no_artificial,
                                         options=['None'] + self.dependent_no_artificial,
                                         value='None',
                                         description='Reference Time Course',
                                         visible=False)
        target_phase_shift_course = Dropdown(values=['None'] + self.dependent_no_artificial,
                                             options=['None'] + self.dependent_no_artificial,
                                             value='None',
                                             description='Calculate Phase Shift for Time Course',
                                             visible=False)

        self.phase_shift_box2.children = [calculate_phase_shift_titration,
                                          reference_time_course,
                                          target_phase_shift_course]

        self.phase_shift_box2.calculate = calculate_phase_shift_titration
        self.phase_shift_box2.target = target_phase_shift_course
        self.phase_shift_box2.reference = reference_time_course

        self.response_times_box.calculate_response_time_titration = calculate_response_time_titration
        self.response_times_box.method = response_time_method_titration
        self.response_times_box.reference = reference_state
        self.response_times_box.ignore_if_less = response_time_ignore_if_less_titration

        if controller.ds.number_conservations != 0:
            self.advanced_titrations.children = [initial_string] + \
                                                [value for key, value in self.initial_concentration2.iteritems()] + \
                                                [self.step] + \
                                                [self.backwards] + \
                                                [self.s_system_steady_states_box] + \
                                                [self.checkbox_custom_points_titration] + \
                                                [self.custom_points_titration] + \
                                                [self.response_times_box] + \
                                                [self.phase_shift_box2] + \
                                                [export_excel_data] + \
                                                [self.titration_function] + \
                                                [self.automatic_scale_button] + \
                                                [self.coordinates] + \
                                                [self.solver]
        else:
            self.advanced_titrations.children = [initial_string] + \
                                                [value for key, value in self.initial_concentration2.iteritems()] + \
                                                [self.step] + \
                                                [self.backwards] + \
                                                [self.s_system_steady_states_box] + \
                                                [self.checkbox_custom_points_titration] + \
                                                [self.custom_points_titration] + \
                                                [self.response_times_box] + \
                                                [self.phase_shift_box2] + \
                                                [export_excel_data] + \
                                                [self.titration_function] + \
                                                [self.automatic_scale_button] + \
                                                [self.coordinates]

        # Trajectories
        self.advanced_trajectories = Box(visible=False)
        self.number_points = Slider(description='# Sampled points over one variable', value=1, min=0, max=20, step=1)
        self.trajectory_from = FloatText(description='Concentration Range from:', value=0.001)
        self.trajectory_to = FloatText(description='Concentration Range to:', value=1)
        self.add_custom = Button(description='Add Custom Point')
        self.include_ss_approximation = Checkbox(description='Show S-Systems Fixed Points', value=False)
        self.include_final_point_marker = Checkbox(description='Show Final Point of Trajectory', value=False)
        self.delete = Button(description='Delete', visible=False)
        self.delete.on_click(self.delete_point)

        self.customized_points_container = HBox()
        self.customized_points_results = Box()
        self.customized_points_copy_results = Box()
        self.customized_points_box = Box()
        self.customized_points_container.children = [self.customized_points_box,
                                                     self.customized_points_results,
                                                     self.customized_points_copy_results,
                                                     ]
        self.customized_points_list = []

        if controller.ds.number_conservations != 0:
            self.advanced_trajectories.children = [self.number_points,
                                                   self.trajectory_from,
                                                   self.trajectory_to,
                                                   self.add_custom,
                                                   self.delete,
                                                   self.customized_points_container,
                                                   self.include_ss_approximation,
                                                   self.include_final_point_marker,
                                                   self.coordinates,
                                                   self.automatic_scale_button,
                                                   self.solver,
                                                   self.dynamical_phenotypes_box
                                                   ]
        else:
            self.advanced_trajectories.children = [self.number_points,
                                                   self.trajectory_from,
                                                   self.trajectory_to,
                                                   self.add_custom,
                                                   self.delete,
                                                   self.customized_points_container,
                                                   self.include_ss_approximation,
                                                   self.include_final_point_marker,
                                                   self.coordinates,
                                                   self.automatic_scale_button,
                                                   self.dynamical_phenotypes_box
                                                   ]

        self.add_custom.on_click(self.add_point)



        ## Event Integration
        self.advanced_event_integration = Box(visible=False)
        self.event_integration_copy_final_box = Box()

        solver = [self.solver] if controller.ds.number_conservations != 0 else []
        n_cycles = FloatText(description='Nr. of Integration Cycles')
        event_type = Dropdown(description='Event Type After Each Cycle',
                              values=['None', 'Modify Initial Conditions'],
                              options=['None', 'Modify Initial Conditions'],
                              value='None')

        modification_type = Dropdown(description='Modification',
                                     values=['None', 'Drop if less than'],
                                     options=['None', 'Drop if less than'],
                                     value='None',
                                     visible=False)

        bounds_dic = dict((i, FloatText(description=i, value=0.001)) for i in self.dependent_no_aux)
        bounds_legend = HTML(value='<b>Bounds for Initial Conditions</b>')

        bounds_box = Box()
        bounds_box.children = [bounds_legend] +\
                              [value for key, value in bounds_dic.iteritems()]
        bounds_box.visible = False

        self.advanced_event_integration.children = [HTML(value='<b>Initial Concentrations</b>')] +\
                                                   [value for key, value in self.initial_concentration.iteritems()] + \
                                                   [self.event_integration_copy_final_box] + \
                                                   [self.coordinates] + \
                                                   [n_cycles] + \
                                                   [event_type] + \
                                                   [modification_type] + \
                                                   [bounds_box] + \
                                                   solver

        self.advanced_event_integration.event_type = event_type
        self.advanced_event_integration.modification_type = modification_type
        self.advanced_event_integration.n_cycles = n_cycles
        self.advanced_event_integration.bounds = bounds_dic
        self.advanced_event_integration.bounds_legend = bounds_legend
        self.advanced_event_integration.bounds_box = bounds_box

        advanced_box.children = [self.advanced_time_courses,
                                 self.advanced_titrations,
                                 self.advanced_trajectories,
                                 self.advanced_event_integration]

        event_type.on_trait_change(self.update_event_box, 'value')
        modification_type.on_trait_change(self.update_event_box, 'value')

        # Add on_trait_change:
        for i in [x_axis, y_axis, plot_type, y_axis_secondary]:
            i.on_trait_change(self.update_field, 'value')
        y_axis.on_trait_change(self.titration_custom_function_update_fields, 'value')
        return ('Full System', full_widget)

    def update_event_box(self, name, value):
        tag1 = self.advanced_event_integration.event_type.value
        tag2 = self.advanced_event_integration.modification_type.value
        if tag1 == 'None' or tag2 == 'None':
            self.advanced_event_integration.modification_type.visible = False
            self.advanced_event_integration.bounds_box.visible = False
        if tag1 == 'Modify Initial Conditions':
            self.advanced_event_integration.modification_type.visible = True
            if tag2 == 'Drop if less than':
                self.advanced_event_integration.bounds_box.visible = True
                self.advanced_event_integration.bounds_legend.value = '<b>Set pools to zero if less than: </b>'

    def update_response_time_box(self, name, value):

        self.response_times_box.method.visible = value
        self.response_times_box.reference.visible = value
        self.response_times_box.ignore_if_less.visible = value and self.response_times_box.method.value == 'custom'

    def update_phase_shift_box(self, name, value):

        self.phase_shift_box2.reference.visible = value
        self.phase_shift_box2.target.visible = value

    def populate_ss_functions(self, b):
        members = self.advanced_time_courses.ss_functions.members
        ss_functions = self.advanced_time_courses.ss_functions

        # create new box with elements of steady state function
        name = Text(description="Name", value='')
        f = Text(description="Function", value='')

        col = Dropdown(description='Color',
                       values=['black', 'red', 'darkgreen', 'darkcyan', 'gold', 'grey'],
                       options=['black', 'red', 'darkgreen', 'darkcyan', 'gold', 'grey'],
                       value='black')

        axis = Dropdown(description='Axis',
                        values=['Primary', 'Secondary'],
                        options=['Primary', 'Secondary'],
                        value='Primary')
        member = Box(children=[name, f, axis, col])
        member.name_ = name
        member.f_ = f
        member.axis_ = axis
        member.col_ = col
        members += [member]

        self.advanced_time_courses.ss_functions.members = members

        delete_button = Button(description="Delete")
        delete_button.on_click(self.decrease_ss_functions)

        ss_functions.children = members + [delete_button]

    def update_visibility(self, name, value, widget=None):
        if widget is not None:
            widget.visible = value
            self.step.visible = not value
            self.x_min.visible = not value
            self.x_max.visible = not value

    def update_visibility2(self, name, value, widget=None):
        if widget is not None:
            widget.visible = value

    def decrease_ss_functions(self, b):
        ss_functions = self.advanced_time_courses.ss_functions
        members = self.advanced_time_courses.ss_functions.members
        delete_button = Button(description="Delete")
        delete_button.on_click(self.decrease_ss_functions)

        if len(members) > 0:
            members = members[:-1]
            self.advanced_time_courses.ss_functions.members = members
            if len(members) > 0:
                ss_functions.children = members + [delete_button]
            else:
                ss_functions.children = members

    def get_multipliers(self, s, s_mult):

        # print("get_multipliers function got s = {} and s_mult = {}".format(s, s_mult))

        output_list = [li for li in difflib.ndiff(s, s_mult) if li[0] != ' ']
        processed_output_list = []
        for element in output_list:
            processed_output_list.append(element.replace('+ ', ''))
        multiplier = "".join(processed_output_list)

        ## delete a * sign at the beginning
        if multiplier[0] == '*':
            multiplier = multiplier[1:]

        # delete * sign in the last position.
        if multiplier[-1] == '*':
            multiplier = multiplier[:-1]

        # replace double ** by single
        multiplier = multiplier.replace('**', '*')

        # split multiplier in multiple elements
        multiplier = multiplier.split('*')

        # loop over each element and replace ^ by ^- and add ^-1 if ^ not present:
        multiplier_mod = []
        for l in multiplier:
            if l.find('^') != -1:
                # replace ^ sign by ^-
                multiplier_mod.append(l.replace('^', '^-'))
            else:
                multiplier_mod.append(l + '^-1')

        multiplier = '*'.join(multiplier_mod)

        # delete * sign in the last position.
        if multiplier[-1] == '*':
            multiplier = multiplier[:-1]
        # add * sign in the first position
        if multiplier[0] != '*':
            multiplier = '*' + multiplier
        return multiplier


    def process_conservation_equations(self, variables):
        # variables is a list of string with variables matching self.conservations.
        self.conservations_dictionary = {}
        self.conservations_equal_zero_dictionary = {}
        self.conservations_variable_list = variables

        if len(variables) != len(self.conservations):
            print("Warning! Conservation relationships won't be correctly parsed. Don't trust results.")

        # print("The variables is ", variables)

        for indx, eq in enumerate(self.conservations):
            parts = str(eq).split("-")

            # print("the parts variable is ", parts)

            # 1. We first verify that the variables[index] is present in the respective conservation and extract the
            # expression in which the variables[index] is present. This is stored in element_to_replace
            element_to_replace = ''
            for element in parts:
                if variables[indx] in dspace.Expression(element).variables:
                    parts_to_replace = variables[indx]
                    element_to_replace = element

            # print("The variable[index] is ", variables[indx])
            # print("The element_to_replace is ", element_to_replace)

            # 2. We check if element_to_replace is different from variables[index]. If so, we need to extract the multipliers
            if element_to_replace == '':
                print("Warning! Conservation relationships won't be correctly parsed. Don't trust results.")
            else:
                if element_to_replace == variables[indx]:
                    multiplier = ''
                else:
                    multiplier = self.get_multipliers(variables[indx], element_to_replace)

            # print("The multipliers are: ", multiplier)

            # 3. Loop over elements of parts, exclude element_to_replace and add multiplier.
            parts_modified = []
            for element in parts:
                if element != element_to_replace:
                        parts_modified.append(element + multiplier)

            # print("parts_modified is: ", parts_modified)
            # print("parts_to_replace is: ", parts_to_replace)

            # 3. We loop over parts_modified and construct the subset parts_filtered, which excludes the element_to_replace
            parts_filtered =[]
            for element in parts_modified:
                parts_filtered.append(element.replace("0", parts_to_replace))

            # print("parts_filtered is ", parts_filtered)

            if len(parts_to_replace) == 0:
                print("the len of parts_to_replace is ", len(parts_to_replace))
                print("Warning, parsing of algebraic constraints is wrong! Do not trust results.")
                print("parts_to_replace  is ", parts_to_replace)
                print("variables are ", variables)

            new_eq = "-".join(parts_filtered)
            # print("The new equation is: ", new_eq)

            self.conservations[indx] = Expression(new_eq)
            self.conservations_dictionary.update({variables[indx]: Expression(new_eq)})
            self.conservations_equal_zero_dictionary.update({variables[indx]: eq})

            # Now we construct a factor dicitonary to take care of coefficients in the LHS
            # Note that the dictionary self.coef_dic is not required anymore.

            coef_list = parts_to_replace[0].split(variables[indx])
            try:
                coef = float(coef_list[0])
            except:
                coef = 1
            self.coef_dic.update({variables[indx]: coef})

    def integrate_conserved_relationships(self, depend_expression):

        delete_list = []
        modified_depend_expression = []
        for conservation in self.conservations:

            # first get the list of dependent variables involved in the conservation
            involved_conserved_variables = [var for var in conservation.variables if var in self.dependent_no_aux]
            unique_involved_conserved_variables = [var for var in involved_conserved_variables if
                                                   var not in delete_list]
            if len(unique_involved_conserved_variables) != 0:
                var_to_delete = unique_involved_conserved_variables[0]
                delete_list.append(var_to_delete)
            else:
                print("Please Check Syntax of Conservation Constraints")

        # Now loop over depend_expression. Delete if left hand side is equal to var_to_delete
        for equation in depend_expression:
            if str(equation.lhs).replace('.', '') not in delete_list:
                modified_depend_expression.append(equation)
        self.process_conservation_equations(delete_list)
        return modified_depend_expression

    def generate_ode_equations(self):

        controller = self.controller

        # construct a list  of len(dependent) first elements of ds.equations.
        depend_expression = controller.ds.equations[0:len(self.dependent_no_aux)]
        number_conservations = controller.ds.number_conservations
        n = len(self.dependent_no_aux)
        aux_expression = controller.ds.equations[n:
                                                 n + len(controller.ds.auxiliary_variables) - number_conservations]

        # use the subst method of the expression object to substitute auxiliary expressions into dependent
        # expressions. first, generate a dictionary out of aux_expression list.
        aux_dic = dict((str(element.lhs), element.rhs) for element in aux_expression)

        # use the subst method on the dictionary itself to take care of nested auxiliary variables present in metabolic models.
        for element in aux_dic:
            aux_dic[element] = str(aux_dic[element].subst(aux_dic))
        self.aux_dic = aux_dic

        # substitute and create list.
        if len(controller.ds.auxiliary_variables) - number_conservations != 0:
            depend_expression = [element.subst(aux_dic) for element in depend_expression]

        # check if there are conservation constraints present. If so, modify the array for depend_expression
        if number_conservations != 0:
            self.conservations = controller.ds.equations[-number_conservations:]
            depend_expression = self.integrate_conserved_relationships(depend_expression)
            self.controller.dae_equations = depend_expression + self.conservations
        self.ode_equations = depend_expression

        # create parameters pool & make it available for all methods.
        parameters = dspace.VariablePool(names=self.dependent_no_aux + controller.ds.independent_variables)
        parameters.update(controller.pvals)
        self.parameters = parameters

        # print("all equations are: ")
        # print(controller.ds.equations)
        # print("depend_expression is:")
        # print(depend_expression)
        # print("The number of dependent_no_aux is ", len(self.dependent_no_aux))
        # print("The number of auxiliary no conserved is ",len(controller.ds.auxiliary_variables) - number_conservations)
        # print("Equations for the auxiliary variables are:")
        # print(aux_expression)
        # print("auxiliary variables are:")
        # print(controller.ds.auxiliary_variables)
        # print("is the system conserved? ", controller.ds.is_conserved)
        # print("the number of conservations is: ", controller.ds.number_conservations)
        # print("depend expression after .subst method:")
        # print(depend_expression)
        # print("The variables contained in the firxt conservation constraint are ", self.conservations[0].variables)
        # print("The ODE equations are: ", self.ode_equations)
        # print("The conserved equations after processing are ", self.conservations)
        # print("The dictionary of coefficients is ", self.coef_dic)
        # print("The conserved equations in zero form are ", self.conservations_equal_zero_dictionary)

    def update_initial_conditions(self, b):

        if self.initial_conditions is not None:
            for key in self.initial_conditions.keys():
                if key in self.dependent_no_aux:
                    self.initial_concentration2[key].value = self.initial_conditions[key]
                    self.initial_concentration[key].value = self.initial_conditions[key]
                else:
                    self.initial_concentration2[key].value = 0.001
                    self.initial_concentration[key].value = 0.001

    def update_scales(self, b):

        if b.description == 'Set Scales':
            b.description = 'Automatic Scales'
            self.automatic_scale = False
            if self.type.value == 'Time Courses':
                self.x_min.visible = False
                self.x_max.visible = False

                self.y_min.description = 'Y-Min'
                # self.y_min.value = 0.01

                self.y_max.description = 'Y-Max'
                # self.y_max.value = 100

                self.y_max.visible = True
                self.y_min.visible = True

                if self.ylabel_secondary.value != 'None' or len(self.advanced_time_courses.ss_functions.members) > 0:
                    self.y_max_secondary.visible = True
                    self.y_min_secondary.visible = True


            elif self.type.value == 'Titrations':
                self.x_min.visible = True
                self.x_min.description = 'X-Min'
                # self.x_min.value = 0.01

                self.x_max.visible = True
                self.x_max.description = 'X-Max'
                # self.x_max.value = 100

                self.y_max.visible = True
                self.y_max.description = 'Y-Max'
                # self.y_max.value = 100

                self.y_min.visible = True
                self.y_min.description = 'Y-Min'
                # self.y_min.value = 0.01
            else:
                self.x_min.visible = True
                self.x_min.description = self.process_label_string('X-Min')
                # self.x_min.value = self.process_results_vector(0.01)

                self.x_max.visible = True
                self.x_max.description = self.process_label_string('X-Max')
                # self.x_max.value = self.process_results_vector(100)

                self.y_max.visible = True
                self.y_max.description = self.process_label_string('Y-Max')
                # self.y_max.value = self.process_results_vector(100)

                self.y_min.visible = True
                self.y_min.description = self.process_label_string('Y-Min')
                # self.y_min.value = self.process_results_vector(0.01)

        else:
            b.description = 'Set Scales'
            self.automatic_scale = True
            if self.type.value == 'Titrations':
                self.x_min.visible = True
                self.x_max.visible = True
                self.y_max.visible = False
                self.y_min.visible = False
            else:
                self.x_min.visible = False
                self.x_max.visible = False
                self.y_max.visible = False
                self.y_min.visible = False
                self.y_max_secondary.visible = False
                self.y_min_secondary.visible = False

    def update_field(self, name, value):

        controller = self.controller

        if 'Custom Function' in self.ylabel.values.keys() and str(self.type.value) != 'Titrations':
            self.ylabel.options = self.dependent_no_artificial
            self.ylabel.values = dict((i, i) for i in self.ylabel.options)

        if value == 'Time Courses':
            self.xlabel.values = {'Time': 'Time'}
            self.title.value = 'Analysis of the ' + controller.name + ' by system design space showing a time course plot'
            if self.ylabel_secondary.value == 'None':
                self.caption.value = 'Concentration of ' + self.ylabel.value + ' is plotted against time.'
            else:
                self.caption.value = 'Concentrations of ' + self.ylabel.value + ' and ' + self.ylabel_secondary.value + ' are plotted against time.'
            # self.x_min.value = 0
            # self.x_max.value = 100
            self.x_min.visible = False
            self.x_max.visible = False
            if hasattr(self, 'y_min'):
                self.y_min.visible = False
                self.y_max.visible = False

            self.advanced_time_courses.visible = True
            self.advanced_titrations.visible = False
            self.advanced_event_integration.visible = False
            self.advanced_trajectories.visible = False
            self.ylabel_secondary.visible = True

        elif value == 'Titrations':
            self.xlabel.values = dict((k, k) for k in controller.ds.independent_variables)
            self.title.value = 'Analysis of the ' + controller.name + ' by system design space showing a titration plot'
            dep_variable = self.ylabel.value if str(self.ylabel.value) != 'Custom Function' else 'Custom Function'
            self.caption.value = 'Concentration of ' + dep_variable + ' is plotted as a function of ' \
                                 + self.xlabel.value + '.'
            self.step.description = 'Number of Points for ' + str(self.xlabel.value)
            self.advanced_time_courses.visible = False
            self.advanced_titrations.visible = True  # if advanced options need to be set, change to True
            self.advanced_event_integration.visible = False
            self.advanced_trajectories.visible = False
            self.x_min.visible = True
            self.x_min.description = 'X-Min'
            # self.x_min.value = 0.01
            self.x_max.visible = True
            self.x_max.description = 'X-Max'
            # self.x_max.value = 100
            # self.y_min.visible = False
            # self.y_max.visible = False
            self.y_max_secondary.visible = False
            self.y_min_secondary.visible = False
            self.ylabel_secondary.visible = False
            self.ylabel.options = self.dependent_no_artificial + ['Custom Function']
            self.ylabel.values = dict((i, i) for i in self.ylabel.options)


        elif value == 'Trajectories':
            self.xlabel.values = dict((k, k) for k in self.dependent_no_aux)
            self.title.value = 'Analysis of the ' + controller.name + ' by system design space showing a trajectories plot'
            self.caption.value = 'Concentration of ' + self.ylabel.value + ' as a function of ' + self.xlabel.value + '.'
            # self.x_min.value = 0.001

            self.advanced_time_courses.visible = False
            self.advanced_titrations.visible = False
            self.advanced_event_integration.visible = False
            self.advanced_trajectories.visible = True
            self.x_min.visible = False
            self.x_max.visible = False
            self.y_min.visible = False
            self.y_min.visible = False
            self.automatic_scale = True
            self.y_max_secondary.visible = False
            self.y_min_secondary.visible = False
            self.ylabel_secondary.visible = False

            self.automatic_scale_button.description = 'Set Scales'

        elif value == 'Event Integration':

            self.xlabel.values = {'Time': 'Time'}
            self.title.value = 'Analysis of the ' + controller.name + ' by system design space showing a time course plot'
            if self.ylabel_secondary.value == 'None':
                self.caption.value = 'Concentration of ' + self.ylabel.value + ' is plotted against time.'
            else:
                self.caption.value = 'Concentrations of ' + self.ylabel.value + ' and ' + self.ylabel_secondary.value + ' are plotted against time.'
            self.x_min.visible = False
            self.x_max.visible = False
            if hasattr(self, 'y_min'):
                self.y_min.visible = False
                self.y_max.visible = False

            self.advanced_event_integration.visible = True
            self.advanced_time_courses.visible = False
            self.advanced_titrations.visible = False
            self.advanced_trajectories.visible = False

        else:
            self.update_field(None, self.type.value)

    def titration_custom_function_update_fields(self, name, value):
        if str(self.ylabel.value) != 'Custom Function':
            self.titration_function.visible = False
            self.coordinates.visible = True
            return
        self.coordinates.visible = False
        self.titration_function.visible = True


    def add_point(self, b):
        controller = self.controller
        # Add element to list self.customized_points_list.
        self.delete.visible = True
        # 1. Generate dictionary
        custom_point = dict(
            (i, FloatText(description=i, value=0.01)) for i in self.dependent_no_aux)
        # 2. Append dictionary to the list: self.customized_points_list
        self.customized_points_list.append(custom_point)
        # 3. refresh widget. Call function self.actualize_custom_box
        self.actualize_custom_box()

    def delete_point(self,b):
        del self.customized_points_list[-1]
        if len(self.customized_points_list) == 0:
            self.delete.visible = False
        else:
            self.delete.visible = True
        self.actualize_custom_box()

    def actualize_custom_box(self):
        # Add children to self.customized_points_box. All customized points are located in a list containing
        # dictionaries. Each element of the list represents a customized point.
        children = []
        # i = 1
        for element in self.customized_points_list:
            # children.append(HTML(value='<b> Point ' + str(i) + '<b>'))
            children += [value for key, value in element.iteritems()]
            # i = i + 1
        self.customized_points_box.children = children

    def call_ode_case(self, b):
        b.pvals = self.controller.pvals.copy()
        if self.type.value == 'Time Courses':
            self.time_course_case(b)
        if self.type.value == 'Titrations':
            self.titrations_case(b)
        if self.type.value == 'Trajectories':
            self.trajectories_case(b)
        if self.type.value == 'Event Integration':
            self.event_integration_case(b)

    def complement_results_dic_with_aux_variables(self, results_dic, n):

        # This function should add vectors for auxiliary variables (algebraic constraints) to the results_dictionary.
        # Expressions defining these variables are defined in self.aux_dic
        aux_dic = self.aux_dic
        # 1. update results_dic with keys of aux_dic:
        keys_aux_dic = aux_dic.keys()
        keys_dep_var = results_dic.keys()
        results_dic.update({key: [] for key in keys_aux_dic})
        p_values = self.controller.pvals.copy()

        # 2. loop over vectors, update dictionary, compute value of auxiliary variables and update corresponding vectors

        for ndx in range(n):
            # update parameter dictionary with values for dependent variables
            p_values.update({key: results_dic[key][ndx] for key in keys_dep_var})
            # compute value of auxiliary variables and update corresponding vectors
            for element in keys_aux_dic:
                results_dic[element].append(Expression(aux_dic[element]).eval_with_values(p_vals=p_values))
        # 3. Transform list into np.arrays for auxiliary variables. This is done for data consistency.
        for element in keys_aux_dic:
            results_dic[element] = np.array(results_dic[element])

    def complement_results_dic(self, results_dic):

        conservation_dic = self.conservations_dictionary
        nr_rows = len(results_dic[results_dic.keys()[0]])
        nr_colums = len(conservation_dic.keys())
        new_results = np.zeros([nr_rows, nr_colums])
        parameters_dic = self.controller.pvals.copy()

        for row, _ in enumerate(results_dic[results_dic.keys()[0]]):
            #update dictionary of parameters
            sub_dic = dict((var, results_dic[var][row]) for var in results_dic.keys())
            parameters_dic.update(sub_dic)
            for col, key in enumerate(conservation_dic.keys()):
                expression = conservation_dic[key]
                new_results[row,col] = expression.rhs.eval_with_values(p_vals=parameters_dic) / self.coef_dic[key]

        # now update results_dic with rows from new_results matrix

        for col, key in enumerate(conservation_dic.keys()):
            results_dic.update({key: new_results[:,col]})

    def calculate_ss_function(self, results_dic):

        # This function should update results_dic with vectors for each member of
        # self.advanced_time_courses.ss_functions.members
        members = self.advanced_time_courses.ss_functions.members
        pvals = self.controller.pvals.copy()

        ll = results_dic[results_dic.keys()[0]]
        for member in members:

            v = []
            # Generate Expression
            if member.f_.value != '' and member.name_.value != '':
                expr = Expression(str(member.f_.value))
            else:
                continue
            for i, _ in enumerate(ll):
                # update dictionary
                sub_dic = dict((var, results_dic[var][i]) for var in results_dic.keys())
                pvals.update(sub_dic)
                # evaluate expression and append.
                v.append(expr.eval_with_values(pvals))
            results_dic.update({str(member.name_.value): v})


    def solve_ode(self, yinit, yinit_conserved, time):

        # This function is used to wrap the solvers ODEint and DAE solvers. It takes as input yinit and time and should
        # return a dictionary containing the solution.
        # important variables are:

        # self.ode_equations --> list of DSExpressions containing differential equations
        # self.conservations_equal_zero_dictionary --> dictionary of expressions containing conservation relationships
        # self.controller.ds.number_conservations --> contain the number of conservations. Equal to 0 if system is not
        # conserved
        # self.solver.value --> contains the solver that should be used to generate the vector y. Possibilities include:
        # options=['IDA', 'RADAU5', 'GLIMDA', 'ODEINT'],
        # self.conservations_variable_list --> list containing strings with conserved variables

        ode = self.ode_equations
        n_conservations = self.controller.ds.number_conservations
        if n_conservations != 0:
            conservations = self.conservations_equal_zero_dictionary
            conservations_list = self.conservations_variable_list
        solver = self.solver.value

        # Standard ODE, no conservations. For titrations, this function returns vector y.
        if n_conservations == 0:
            # rtol = 1e-10  #defaults are 1.49012e-8
            # atol = 1e-10  #defaults are 1.49012e-8
            y = odeint(self.construct_fun, yinit, time)
            if self.type.value == 'Titrations':
                return y
            results_dic = {}
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                results_dic.update({var: y[:, counter]})
            # if len(self.advanced_time_courses.ss_functions.members) > 0:
            #     self.calculate_ss_function(results_dic)
            return results_dic

        # Conservations using standard solver ODEINT. This option is not recommended.
        elif n_conservations != 0 and solver == 'ODEINT':
            y = odeint(self.construct_fun, yinit, time)
            if self.type.value == 'Titrations':
                return y
            results_dic = {}
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                results_dic.update({var: y[:, counter]})
            if self.conservations_dictionary is not None:
                self.complement_results_dic(results_dic)
            return results_dic

        # Recommended option to solve DAE systems.
        elif n_conservations != 0 and solver != 'ODEINT':

            # 0. Expand y0 to include algebraic constraints
            y0 = np.append(yinit, yinit_conserved)

            # 1. calculate yd0 from yd
            yd0 = self.calculate_yd0_from_y0(yinit, yinit_conserved)

            # 2. create implicit problem
            imp_mod = Implicit_Problem(self.construct_fun_DAE, y0, yd0)

            # 3. set the algebraic components
            v = np.zeros(len(ode) + n_conservations) + 1
            v[-n_conservations:] = 0
            imp_mod.algvar = v

            if solver == 'IDA':
                # Create an Assimulo implicit solver (IDA)
                imp_sim = IDA(imp_mod)  # Create a IDA solver

                # Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
                imp_sim.make_consistent('IDA_YA_YDP_INIT')

            elif solver == 'RADAU5':
                imp_sim = Radau5DAE(imp_mod)

            # Sets the paramters
            imp_sim.atol = 1e-6  # Default 1e-6
            imp_sim.rtol = 1e-6  # Default 1e-6
            imp_sim.suppress_alg = False  # Suppress the algebraic variables on the error test
            imp_sim.verbosity = 50 # Default 50, scream = 10
            imp_sim.maxsteps = 10000  # default 10000

            # Simulate
            t, y, yd = imp_sim.simulate(time[-1], len(time)-1)
            if self.type.value == 'Titrations':
                return y

            # pack results in a dictionary. First get vectors for ode_variables
            results_dic = {}
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                results_dic.update({var: y[:, counter]})
            # now get vectors for algebraic constraints
            for counter2, element in enumerate(conservations_list):
                n = counter + counter2 + 1
                results_dic.update({element: y[:, n]})
            # if len(self.advanced_time_courses.ss_functions.members) > 0:
            #     self.calculate_ss_function(results_dic)
            return results_dic

    def solve_ode_cycles(self, yinit, yinit_conserved, time, n_cycles, bounds, comparison_type):

        # This function is used to wrap the solvers ODEint and DAE solvers. It takes as input yinit and time and should
        # return a dictionary containing the solution. The behavior of the integration is defined by the variables
        # n_cycles, bounds, comparison_type

        # important variables are:

        # self.ode_equations --> list of DSExpressions containing differential equations
        # self.conservations_equal_zero_dictionary --> dictionary of expressions containing conservation relationships
        # self.controller.ds.number_conservations --> contain the number of conservations. Equal to 0 if system is not
        # conserved
        # self.solver.value --> contains the solver that should be used to generate the vector y. Possibilities include:
        # options=['IDA', 'RADAU5', 'GLIMDA', 'ODEINT'],
        # self.conservations_variable_list --> list containing strings with conserved variables

        ode = self.ode_equations
        n_conservations = self.controller.ds.number_conservations
        if n_conservations != 0:
            conservations = self.conservations_equal_zero_dictionary
            conservations_list = self.conservations_variable_list
        solver = self.solver.value

        # Standard ODE, no conservations. For titrations, this function returns vector y.
        if n_conservations == 0:

            # we need to run odeint for a number of times defined by n_cycles
            for n in range(int(n_cycles)):

                y = odeint(self.construct_fun, yinit, time)
                # update yinit as final values of y
                for counter, element in enumerate(self.ode_equations):
                    yinit[counter] = y[-1, counter]
                # Compare yinit with bounds and drop if comparison_type is less than
                sum_yinit = 0
                if comparison_type == 'Drop if less than':
                    for i in range(len(self.ode_equations)):
                        yinit[i] = yinit[i] if yinit[i] > bounds[i] else 0
                        sum_yinit += yinit[i]
                    # rescale yinit
                    for i in range(len(self.ode_equations)):
                        yinit[i] = yinit[i]/sum_yinit

            if self.type.value == 'Titrations':
                return y
            results_dic = {}
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                results_dic.update({var: y[:, counter]})

            return results_dic

        # Recommended option to solve DAE systems.
        elif n_conservations != 0 and solver != 'ODEINT':
            print("Currently not supported")

            # # 0. Expand y0 to include algebraic constraints
            # y0 = np.append(yinit, yinit_conserved)
            #
            # # 1. calculate yd0 from yd
            # yd0 = self.calculate_yd0_from_y0(yinit, yinit_conserved)
            #
            # # 2. create implicit problem
            # imp_mod = Implicit_Problem(self.construct_fun_DAE, y0, yd0)
            #
            # # 3. set the algebraic components
            # v = np.zeros(len(ode) + n_conservations) + 1
            # v[-n_conservations:] = 0
            # imp_mod.algvar = v
            #
            # if solver == 'IDA':
            #     # Create an Assimulo implicit solver (IDA)
            #     imp_sim = IDA(imp_mod)  # Create a IDA solver
            #
            #     # Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
            #     imp_sim.make_consistent('IDA_YA_YDP_INIT')
            #
            # elif solver == 'RADAU5':
            #     imp_sim = Radau5DAE(imp_mod)
            #
            # # Sets the paramters
            # imp_sim.atol = 1e-6  # Default 1e-6
            # imp_sim.rtol = 1e-6  # Default 1e-6
            # imp_sim.suppress_alg = False  # Suppress the algebraic variables on the error test
            # imp_sim.verbosity = 50  # Default 50, scream = 10
            # imp_sim.maxsteps = 10000  # default 10000
            #
            # # Simulate
            # t, y, yd = imp_sim.simulate(time[-1], len(time) - 1)
            # if self.type.value == 'Titrations':
            #     return y
            #
            # # pack results in a dictionary. First get vectors for ode_variables
            # results_dic = {}
            # for counter, element in enumerate(self.ode_equations):
            #     var = str(element.lhs).replace('.', '')
            #     results_dic.update({var: y[:, counter]})
            # # now get vectors for algebraic constraints
            # for counter2, element in enumerate(conservations_list):
            #     n = counter + counter2 + 1
            #     results_dic.update({element: y[:, n]})
            # # if len(self.advanced_time_courses.ss_functions.members) > 0:
            # #     self.calculate_ss_function(results_dic)
            # return results_dic


    def calculate_yd0_from_y0(self, yinit, yinit_conserved):

        yd0 = []

        # create dictionary of parameters and update it by taking values from toolbox
        par = self.parameters

        # update dictionary of parameter values to consider initial conditions for Xd_t pools
        for indx, element in enumerate(self.ode_equations):
            par[str(element.lhs).replace('.', '')] = yinit[indx]

        # update dictionary of parameter values to consider intial conditions for Xd_a_c pools
        for indx, element in enumerate(self.conservations_variable_list):
            par[element] = yinit_conserved[indx]

        # Now start populating yd0. Start evaluating the rhs of self.ode_equations
        for eq in self.ode_equations:
            yd0.append(eq.rhs.eval_with_values(p_vals=par))

        # now append a list of zeros to account for conserved variables
        for eq in self.conservations_variable_list:
            yd0.append(0)

        return yd0

    def time_course_case(self, b):
        # set initial conditions.
        yinit = np.array([])
        yinit_conserved = np.array([])
        # print("The system for the ODE solver is:")

        initial_conditions_string = ''

        for element in self.ode_equations:
            # print(element)
            yinit = np.append(yinit, self.initial_concentration[str(element.lhs).replace('.', '')].value)
            initial_conditions_string += str(element.lhs).replace('.', '') + ' = ' + \
                                         str(self.initial_concentration[str(element.lhs).replace('.', '')].value) + '; '

        if self.controller.ds.number_conservations != 0:
            for element in self.conservations_variable_list:
                yinit_conserved = np.append(yinit_conserved, self.initial_concentration[element].value)

        time = np.linspace(self.x_min_time.value, self.x_max_time.value, self.step_time.value)

        # attach string containing initial conditions to self.controller object
        self.controller.initial_conditions_caption = initial_conditions_string[:-2]

        # Actualize parameters
        self.parameters.update(self.controller.pvals)
        results_dic = self.solve_ode(yinit, yinit_conserved, time)
        # complement results_dic with values for auxiliary variables (algebraic constraints)
        self.complement_results_dic_with_aux_variables(results_dic, len(time))
        # compute additional ss functions
        if len(self.advanced_time_courses.ss_functions.members) > 0:
            self.calculate_ss_function(results_dic)

        # print results_dic
        self.make_plot(time, results_dic, b)
        copy_steady_state_button = Button(description='Copy Final Steady States')
        copy_steady_state_button.results = results_dic
        copy_steady_state_button.case = 'Time_Course'
        copy_steady_state_button.on_click(self.copy_steady_states)

        export_steady_state_button = Button(description='Export Results to .xlsx File')
        export_steady_state_button.results = results_dic
        export_steady_state_button.t = time
        export_steady_state_button.case = 'Time_Course'
        export_steady_state_button.on_click(self.export_numerical_results)

        self.time_courses_copy_final_box.children = [copy_steady_state_button,
                                                     export_steady_state_button]

        if self.include_dynamical_in_excel.value is True:

            xi = dspace.VariablePool(names=self.controller.ds.independent_variables)
            xi.update(self.controller.pvals)
            dom_sig = self.controller.ds.dominant_signature(xi, results_dic, len(time))
            ## Now let's construct a vector containing the case number.
            dynamic_case = []
            for element in dom_sig:
                element = element.strip().replace(' ', '').replace(' ', '').replace(' ', '')
                dynamic_case.append(self.controller.ds.case_number_for_signature(element))

            ## Let's attach these variables to the self.dynamical_phenotypes_box
            self.dynamical_phenotypes_box.dom_sig = dom_sig
            self.dynamical_phenotypes_box.dynamic_case = dynamic_case

            # Plot dynamical phenotyopes
            if self.include_dynamical_in_plot.value is True:
                self.make_plot_dynamic_phenotypes(time, results_dic, b)

    def event_integration_case(self, b):

        n_cycles = self.advanced_event_integration.n_cycles.value
        if n_cycles == 0:
            return

        # set initial conditions.
        yinit = np.array([])
        yinit_conserved = np.array([])
        # print("The system for the ODE solver is:")

        initial_conditions_string = ''

        for element in self.ode_equations:
            yinit = np.append(yinit, self.initial_concentration[str(element.lhs).replace('.', '')].value)
            initial_conditions_string += str(element.lhs).replace('.', '') + ' = ' + \
                                         str(self.initial_concentration[str(element.lhs).replace('.', '')].value) + '; '

        if self.controller.ds.number_conservations != 0:
            for element in self.conservations_variable_list:
                yinit_conserved = np.append(yinit_conserved, self.initial_concentration[element].value)

        time = np.linspace(self.x_min_time.value, self.x_max_time.value, self.step_time.value)

        # attach string containing initial conditions to self.controller object
        self.controller.initial_conditions_caption = initial_conditions_string[:-2]

        # Actualize parameters
        self.parameters.update(self.controller.pvals)


        comparison_type = self.advanced_event_integration.modification_type.value

        bounds = []
        bounds_dic = self.advanced_event_integration.bounds

        for counter, element in enumerate(self.ode_equations):
            var = str(element.lhs).replace('.', '')
            bounds.append(bounds_dic[var].value)

        results_dic = self.solve_ode_cycles(yinit, yinit_conserved, time, n_cycles, bounds, comparison_type)

        # complement results_dic with values for auxiliary variables (algebraic constraints)
        self.complement_results_dic_with_aux_variables(results_dic, len(time))
        # # compute additional ss functions
        # if len(self.advanced_time_courses.ss_functions.members) > 0:
        #     self.calculate_ss_function(results_dic)

        # print results_dic
        self.make_plot(time, results_dic, b)
        copy_steady_state_button = Button(description='Copy Final Steady States')
        copy_steady_state_button.results = results_dic
        copy_steady_state_button.case = 'Time_Course'
        copy_steady_state_button.on_click(self.copy_steady_states)

        export_steady_state_button = Button(description='Export Results to .xlsx File')
        export_steady_state_button.results = results_dic
        export_steady_state_button.t = time
        export_steady_state_button.case = 'Time_Course'
        export_steady_state_button.on_click(self.export_numerical_results)

        self.event_integration_copy_final_box.children = [copy_steady_state_button,
                                                          export_steady_state_button]


    def titrations_case(self, b):

        # The variables that are important here are:
        # dictionary self.initial_concentration: contains initial concentrations.
        # self.step contain the number of points in the time vector. We will proceed as follows:
        # calculate steady state concentrations for the given set of initial concentrations and initial value of target
        # independent variable
        # parse initial concentrations:

        controller = self.controller
        initial_conditions_string = ''
        yinit = np.array([])
        yinit_conserved = np.array([])

        # extract variables for the calculation of the response times
        calculate_response_time_w = self.response_times_box.calculate_response_time_titration
        method_response_time = self.response_times_box.method
        reference_response_time = self.response_times_box.reference
        ignore_if_less = self.response_times_box.ignore_if_less

        for element in self.ode_equations:
            yinit = np.append(yinit, self.initial_concentration2[str(element.lhs).replace('.', '')].value)
            initial_conditions_string += str(element.lhs).replace('.', '') + ' = ' + \
                                         str(self.initial_concentration2[str(element.lhs).replace('.', '')].value) + '; '

        if controller.ds.number_conservations != 0:
            for element in self.conservations_variable_list:
                yinit_conserved = np.append(yinit_conserved, self.initial_concentration2[element].value)

        self.controller.initial_conditions_caption = initial_conditions_string[:-2]

        # Get the constant initial point for the calculation of response times for initial x_point method
        if calculate_response_time_w.value is True:
            if reference_response_time.value == 'Initial Point':

                yinit_initial_x_point = np.array([])
                yinit_conserved_initial_x_point = np.array([])

                for element in self.ode_equations:
                    yinit_initial_x_point = np.append(yinit_initial_x_point,
                                                      self.initial_concentration2[
                                                          str(element.lhs).replace('.', '')].value)

                if controller.ds.number_conservations != 0:
                    for element in self.conservations_variable_list:
                        yinit_conserved_initial_x_point = np.append(yinit_conserved_initial_x_point,
                                                                    self.initial_concentration2[element].value)

        # Generate vector containing values for the selected variable on the x-axis.
        target_x_parameter = np.linspace(np.log10(self.x_min.value), np.log10(self.x_max.value), num=self.step.value)
        target_x_parameter = 10 ** target_x_parameter

        # update target_x_parameter as a function of self.checkbox_custom_points_titration
        if self.checkbox_custom_points_titration.value is True:
            aux = self.custom_points_titration.value.split(',')
            aux_val = []
            for element in aux:
                aux_val.append(float(element))
            target_x_parameter = aux_val

        # Actualize parameters
        self.parameters.update(self.controller.pvals)
        self.parameters.update({self.xlabel.value: target_x_parameter[0]})


        # generate time vector.
        time = np.linspace(self.x_min_time.value, self.x_max_time.value, self.step_time.value)

        # Generate first steady state point
        y = self.solve_ode(yinit, yinit_conserved, time)

        # parse solution.
        # create a dictionary with data.
        results_dic = {}
        for counter, element in enumerate(self.ode_equations):
            results_dic.update({str(element.lhs).replace('.', ''): [y[-1, counter]]})

        if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
            for counter2, element in enumerate(self.conservations_variable_list):
                results_dic.update({element: [y[-1, counter+counter2+1]]})

        # Initialize and parse dictionaries for response times.
        if calculate_response_time_w.value is True:

            method = method_response_time.value
            threshold = ignore_if_less.value / 100

            if method_response_time.value == 'custom':
                indx = -1
            else:
                indx = 0

            results_dic_response = {}
            # Generate response times
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                vector = y[:, counter]
                response = calculate_response_time_vector(time, vector, method, threshold)[indx]
                results_dic_response.update({var: [response]})

            if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                for counter2, var in enumerate(self.conservations_variable_list):
                    vector = y[:, counter+counter2+1]
                    response = calculate_response_time_vector(time, vector, method, threshold)[indx]
                    results_dic_response.update({var: [response]})

        # initialize and parse dictionaries for phase shift and amplitude
        el = self.phase_shift_box2
        if el.calculate.value is True and el.target.value != 'None' and el.reference.value != 'None':
            results_dic_PhaseShift = {}
            results_dic_Amplitude = {}
            variable_list = [str(element.lhs).replace('.', '') for element in self.ode_equations]
            reference_vector = find_vector(variable_list, y, el.reference.value)
            target_vector = find_vector(variable_list, y, el.target.value)
            phase_shift, amplitude = calculate_phase_shifts_from_vectors(time,
                                                                         reference_vector,
                                                                         target_vector,
                                                                         mode='Excel')
            results_dic_PhaseShift.update({el.target.value: [phase_shift]})
            results_dic_Amplitude.update({el.target.value: [amplitude]})

        # loop over other elements of target_x_parameter
        for w in target_x_parameter[1:]:
            self.parameters.update({self.xlabel.value: w})
            # initial concentrations correspond to last steady state point in results_dic.
            # Remember that the order of the dictionary does not correspond to the order with wich the elements were stored.
            yinit = np.array([])
            yinit_conserved = np.array([])

            for element in self.ode_equations:
                yinit = np.append(yinit, results_dic[str(element.lhs).replace('.', '')][-1])

            if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                for element in self.conservations_variable_list:
                    yinit_conserved = np.append(yinit_conserved, results_dic[element][-1])

            # Generate next steady state point
            y = self.solve_ode(yinit, yinit_conserved, time)

            # parse and append
            for counter, element in enumerate(self.ode_equations):
                results_dic[str(element.lhs).replace('.', '')].append(y[-1, counter])
            if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                for counter2, element in enumerate(self.conservations_variable_list):
                    results_dic[element].append(y[-1, counter+counter2+1])

            # calculate phase shift and frequencies
            if el.calculate.value is True and el.target.value != 'None' and el.reference.value != 'None':
                reference_vector = find_vector(variable_list, y, el.reference.value)
                target_vector = find_vector(variable_list, y, el.target.value)
                phase_shift, amplitude = calculate_phase_shifts_from_vectors(time,
                                                                             reference_vector,
                                                                             target_vector,
                                                                             mode='Excel')
                results_dic_PhaseShift[el.target.value].append(phase_shift)
                results_dic_Amplitude[el.target.value].append(amplitude)

            # if calculating response times and the reference is initial x_point, update y_init and y_init_conserved
            if calculate_response_time_w.value is True:
                if reference_response_time.value == 'Initial Point':
                    y = self.solve_ode(yinit_initial_x_point, yinit_conserved_initial_x_point, time)

                # parse and append
                method = method_response_time.value
                threshold = ignore_if_less.value / 100

                for counter, element in enumerate(self.ode_equations):
                    var = str(element.lhs).replace('.', '')
                    vector = y[:, counter]
                    response = calculate_response_time_vector(time, vector, method, threshold)[indx]
                    results_dic_response[var].append(response)

                if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                    for counter2, var in enumerate(self.conservations_variable_list):
                        vector = y[:, counter + counter2 + 1]
                        response = calculate_response_time_vector(time, vector, method, threshold)[indx]
                        results_dic_response[var].append(response)

        if self.conservations_dictionary is not None and self.solver.value == 'ODEINT':
            self.complement_results_dic(results_dic)

        # complement dictionary with auxiliary variables
        self.complement_results_dic_with_aux_variables(results_dic, len(target_x_parameter))

        if self.backwards.value is True:
            # now let's go backwards. ##################################################

            # parse initial concentrations:
            yinit = np.array([])
            yinit_conserved = np.array([])
            for element in self.ode_equations:
                yinit = np.append(yinit, results_dic[str(element.lhs).replace('.', '')][-1])
            if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                for element in self.conservations_variable_list:
                    yinit_conserved = np.append(yinit_conserved, results_dic[element][-1])

            # generate vector of parameters
            target_x_parameter_reversed = np.linspace(np.log10(self.x_max.value), np.log10(self.x_min.value), num=self.step.value)
            target_x_parameter_reversed = 10 ** target_x_parameter_reversed

            # update target_x_parameter as a function of self.checkbox_custom_points_titration
            if self.checkbox_custom_points_titration.value is True:
                aux = self.custom_points_titration.value.split(',')
                aux_val = []
                for element in reversed(aux):
                    aux_val.append(float(element))
                target_x_parameter_reversed = aux_val

            # Actualize parameters
            self.parameters.update(self.controller.pvals)
            self.parameters.update({self.xlabel.value: target_x_parameter_reversed[0]})

            # calculate steady state concentration first point.
            y = self.solve_ode(yinit, yinit_conserved, time)
            # y = odeint(self.construct_fun, yinit, time)

            # parse solution.
            # create a dictionary with data.
            results_dic_back = {}
            for counter, element in enumerate(self.ode_equations):
                results_dic_back.update({str(element.lhs).replace('.', ''): [y[-1, counter]]})
            if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                for counter2, element in enumerate(self.conservations_variable_list):
                    results_dic_back.update({element: [y[-1, counter + counter2 + 1]]})

            # loop over other elements.
            for w in target_x_parameter_reversed[1:]:
                self.parameters.update({self.xlabel.value: w})
                # print self.parameters

                # initial concentrations correspond to last steady state point in results_dic
                yinit = np.array([])
                yinit_conserved = np.array([])
                for element in self.ode_equations:
                    yinit = np.append(yinit, results_dic[str(element.lhs).replace('.', '')][-1])

                if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                    for element in self.conservations_variable_list:
                        yinit_conserved = np.append(yinit_conserved, results_dic[element][-1])

                # Generate next steady state point
                y = self.solve_ode(yinit, yinit_conserved, time)
                # y = odeint(self.construct_fun, yinit, time)

                # parse and append
                for counter, element in enumerate(self.ode_equations):
                    results_dic_back[str(element.lhs).replace('.', '')].append(y[-1, counter])
                if self.solver.value != 'ODEINT' and self.controller.ds.number_conservations != 0:
                    for counter2, element in enumerate(self.conservations_variable_list):
                        results_dic_back[element].append(y[-1, counter + counter2 + 1])

            if self.conservations_dictionary is not None and self.solver.value == 'ODEINT':
                self.complement_results_dic(results_dic_back)

            # complement dictionary with auxiliary variables
            self.complement_results_dic_with_aux_variables(results_dic_back, len(target_x_parameter))
            backward = self.process_titration_vector(results_dic_back,
                                                     titration_vector=target_x_parameter_reversed,
                                                     titration_variable=str(self.xlabel.value)
                                                     )

            # End of backwards routine

        # lets generate hin vector
        forward = self.process_titration_vector(results_dic,
                                                titration_vector=target_x_parameter,
                                                titration_variable=str(self.xlabel.value))

        # plotting routine
        fig = plt.figure(figsize=[5.5698, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.1714, 0.2, 0.6, 0.7])

        ax.set_title('Titration Plot')
        ax.set_xlabel('Log('+self.xlabel.value+')')
        ax.plot(np.log10(target_x_parameter), forward, 'blue', linewidth=1)

        if self.backwards.value is True:
            ax.plot(np.log10(target_x_parameter_reversed), backward, 'orange', linewidth=1)

        ax.ticklabel_format(useOffset=False)
        if str(self.ylabel.value) != 'Custom Function':
            fun = self.process_label_string(self.ylabel.value, function=True)
            ax.set_ylabel(self.process_label_string(self.ylabel.value))
        else:
            fun = str(self.titration_function.value)
            ax.set_ylabel(fun)
        rangex = [self.x_min.value, self.x_max.value]
        if self.s_system_steady_states_box.show_ssystems.value is True:
            included_cases_string = str(self.s_system_steady_states_box.only_cases.value)
            if included_cases_string == '':
                included_cases = None
            else:
                included_cases = included_cases_string.strip().replace(' ', '').replace('', '').replace(' ', '').split(',')
            controller.ds.draw_1D_ss_function(ax, fun, controller.pvals,
                                              str(self.xlabel.value),
                                              rangex, ylim=None,
                                              resolution=100,
                                              included_cases=included_cases)
        if self.automatic_scale is False:
            if str(self.ylabel.value) != 'Custom Function':
                ax.set_ylim(self.process_results_vector(self.y_min.value),
                            self.process_results_vector(self.y_max.value))
            else:
                ax.set_ylim(self.y_min.value,
                            self.y_max.value)


        buf = cStringIO.StringIO()
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(buf)
        data = buf.getvalue()
        plt.close()
        # Generate Popup containing figure.
        controller.figures.add_figure(data,
                                      title=self.title.value,
                                      caption=self.caption.value + ' Blue line represents titration for '
                                                                   'increasing values of this variable, while orange '
                                                                   'line represents decreasing values. Black solid line '
                                                                   'is the S-system approximation. ',
                                      pvals=controller.pvals.copy(),
                                      colors=None)

        # Export titrations data to Excel.
        self.advanced_titrations.export_excel_data.visible = True
        self.advanced_titrations.export_excel_data.results_forward = results_dic

        if el.calculate.value is True and el.target.value != 'None' and el.reference.value != 'None':
            self.advanced_titrations.export_excel_data.phase_shift = results_dic_PhaseShift
            self.advanced_titrations.export_excel_data.amplitude = results_dic_Amplitude
        if calculate_response_time_w.value is True:
            self.advanced_titrations.export_excel_data.response_time = results_dic_response
        if self.backwards.value is True:
            self.advanced_titrations.export_excel_data.results_backwards = results_dic_back
            self.advanced_titrations.export_excel_data.independent_parameter_backwards = target_x_parameter_reversed
        self.advanced_titrations.export_excel_data.independent_parameter_forward = target_x_parameter
        self.advanced_titrations.export_excel_data.independent_parameter_identity = self.xlabel.value
        self.advanced_titrations.export_excel_data.case = 'Titrations'
        self.advanced_titrations.export_excel_data.on_click(self.export_numerical_results)

    def process_titration_vector(self, results_dic, titration_vector=None, titration_variable=None):
        fun = self.ylabel.value
        if str(fun) != 'Custom Function':
            vector = results_dic[self.ylabel.value]
            processed_vector = self.process_results_vector(vector)
        else:
            fun = str(self.titration_function.value)
            processed_vector = []
            p_vals = VariablePool()
            p_vals.update(self.controller.pvals)
            if isinstance(fun, Expression):
                expr = fun
            else:
                expr = Expression(fun)
            for n, val in enumerate(results_dic[results_dic.keys()[0]]):
                for key in results_dic.keys():
                    p_vals[key] = results_dic[key][n]
                if titration_vector is not None and titration_variable is not None:
                    p_vals[titration_variable] = titration_vector[n]
                else:
                    raise ValueError('Titration vector or titration variable not provided!')
                value = expr.eval_with_values(p_vals=p_vals)
                processed_vector.append(value)
            results_dic.update({fun: processed_vector})
        return processed_vector

    def trajectories_case(self, b):
        controller = self.controller
        time = np.linspace(self.x_min_time.value, self.x_max_time.value, self.step_time.value)
        fig = plt.figure(figsize=[5.5698, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.1714, 0.2, 0.6, 0.7])
        cap = ""

        # This plot occurs in two steps. First, trajectories for normal distributed points are generated
        # (relevant variables here are self.number_points, self.trajectory_from, self.trajectory_to). Then,
        # tarjectories for single points are generated (self.customized_points_list)

        results_list_sampled_points = []
        results_list = []

        dom_sig_sampled = []
        dynamic_case_sampled = []

        dom_sig_custom = []
        dynamic_case_custom = []

        if self.number_points.value != 0:
            # define vector of initial concentrations.
            initialC = np.linspace(np.log10(self.trajectory_from.value), np.log10(self.trajectory_to.value), num=self.number_points.value)
            initialC = 10**initialC

            # Generate list containing sampled concentrations
            if self.solver.value == 'ODEINT':
                sampled_points = list(itertools.product(initialC, repeat=len(self.dependent_no_aux_no_conserved)))
            else:
                sampled_points = list(itertools.product(initialC, repeat=len(self.dependent_no_aux)))

            #Generate trajectories for each element in sampled_points.

            # Actualize parameters
            self.parameters.update(self.controller.pvals)
            for element in sampled_points:
                yinit = element[:len(self.ode_equations)]
                if controller.ds.number_conservations != 0:
                    yinit_conserved = element[-len(self.conservations_variable_list):]
                else:
                    yinit_conserved = []
                results_dic = self.solve_ode(yinit, yinit_conserved, time)
                self.complement_results_dic_with_aux_variables(results_dic, len(time))
                results_list_sampled_points.append(results_dic)

                ## Get dominat signatures and dynamic phenotypes
                xi = dspace.VariablePool(names=controller.ds.independent_variables)
                xi.update(controller.pvals)
                dom_sig_aux = controller.ds.dominant_signature(xi, results_dic, len(time))
                dom_sig_sampled.append(dom_sig_aux)

                ## Now let's construct a vector containing the case number.
                dynamic_case_aux = []
                for element_ in dom_sig_aux:
                    element_ = element_.strip().replace(' ', '').replace(' ', '').replace(' ', '')
                    dynamic_case_aux.append(self.controller.ds.case_number_for_signature(element_))
                dynamic_case_sampled.append(dynamic_case_aux)

                ax.plot(self.process_results_vector(results_dic[self.xlabel.value]),
                        self.process_results_vector(results_dic[self.ylabel.value]), 'grey', linewidth=1)
                if self.include_final_point_marker.value is True:
                    ax.plot(self.process_results_vector(results_dic[self.xlabel.value])[-1],
                            self.process_results_vector(results_dic[self.ylabel.value])[-1], 'grey',
                            marker="*", markersize=6)

        # Now generate trajectories for each custom point
        if len(self.customized_points_list) != 0:
            self.parameters.update(self.controller.pvals)

            for element in self.customized_points_list:
                # parse initial concentration vector
                yinit = np.array([])
                yinit_conserved = np.array([])

                for variable in self.ode_equations:
                    yinit = np.append(yinit, element[str(variable.lhs).replace('.', '')].value)
                if controller.ds.number_conservations != 0:
                    for variable in self.conservations_variable_list:
                        yinit_conserved = np.append(yinit_conserved, element[variable].value)
                else:
                    yinit_conserved = []
                results_dic = self.solve_ode(yinit, yinit_conserved, time)
                self.complement_results_dic_with_aux_variables(results_dic, len(time))

                ## Get dominat signatures and dynamic phenotypes
                xi = dspace.VariablePool(names=controller.ds.independent_variables)
                xi.update(controller.pvals)
                dom_sig_custom.append(controller.ds.dominant_signature(xi, results_dic, len(time)))

                ## Now let's construct a vector containing the case number.
                dynamic_case_aux = []
                for element_ in dom_sig_sampled:
                    element_ = element_.strip().replace(' ', '').replace(' ', '').replace(' ', '')
                    dynamic_case_aux.append(controller.ds.case_number_for_signature(element_))
                dynamic_case_custom.append(dynamic_case_aux)

                # plot
                results_list.append(results_dic)
                ax.plot(self.process_results_vector(results_dic[self.xlabel.value]),
                        self.process_results_vector(results_dic[self.ylabel.value]), 'blue', linewidth=1)
                if self.include_final_point_marker.value is True:
                    ax.plot(self.process_results_vector(results_dic[self.xlabel.value])[-1],
                            self.process_results_vector(results_dic[self.ylabel.value])[-1],
                            'blue', marker="*", markersize=6)

            self.customized_points_copy_results.results = results_list
            copy_final_steady_states_customized = Button(description='Copy Final Steady States')
            self.customized_points_copy_results.children = [copy_final_steady_states_customized]

            self.customized_points_copy_results.children[0].case = 'Trajectories'
            self.customized_points_copy_results.children[0].on_click(self.copy_steady_states)

        export_steady_state_button = Button(description='Export Results to .xlsx File')
        export_steady_state_button.results = results_list_sampled_points + results_list
        self.dynamical_phenotypes_box.dom_sig = dom_sig_custom + dom_sig_sampled
        self.dynamical_phenotypes_box.dynamic_case = dynamic_case_custom + dynamic_case_sampled
        export_steady_state_button.t = time
        export_steady_state_button.case = 'Trajectories'
        export_steady_state_button.on_click(self.export_numerical_results)

        if len(self.customized_points_list) != 0:
            self.customized_points_copy_results.children = [copy_final_steady_states_customized] + \
                                                       [export_steady_state_button]
        else:
            self.customized_points_copy_results.children = [export_steady_state_button]

        if self.include_ss_approximation.value is True:
            steady_states = self.get_fixed_points_ssystems()
            for s in steady_states:
                x = self.process_results_vector(s[self.xlabel.value])
                y = self.process_results_vector(s[self.ylabel.value])
                col = s["color"]
                ax.plot(x, y, col, marker='o', markersize=6)

        if self.automatic_scale is False:
            ax.set_xlim(self.process_results_vector(self.x_min.value), self.process_results_vector(self.x_max.value))
            ax.set_ylim(self.process_results_vector(self.y_min.value), self.process_results_vector(self.y_max.value))
            #print("setting boundaries")

        ax.set_title('Trajectories Plot')
        ax.set_xlabel(self.process_label_string(self.xlabel.value))
        ax.set_ylabel(self.process_label_string(self.ylabel.value))

        buf = cStringIO.StringIO()
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(buf)
        data = buf.getvalue()
        plt.close()

        if self.include_final_point_marker.value is True:
            cap += " Final point of each trajectory is shown by a star marker. "
        if self.include_ss_approximation.value is True:
            cap += " Steady state of valid phenotypes is shown by a color-coded circle to characterize its stability. " \
                   "Black: a s-system with 0 eigenvalues with (+) real part. Blue: a s-system with 1 eigenvalue with " \
                   "(+) real part. Orange: a s-system with 2 eigenvalues with (+) real part. "

        controller.figures.add_figure(data,
                                      title=self.title.value,
                                      caption=self.caption.value + cap,
                                      pvals=controller.pvals.copy(),
                                      colors=None)

    def get_fixed_points_ssystems(self):

        pvals = self.controller.pvals
        ds = self.controller.ds
        cases = ds.valid_cases(p_bounds=pvals, strict=False)

        ss_steady_states = []

        for i in cases:
            ss_steady_state = ds(i).steady_state(pvals)
            nr_pos_eigen = str(ds(i).positive_roots(pvals)).split("*")[0]
            if nr_pos_eigen == "0":
                col = "black"
            elif nr_pos_eigen == "1":
                col = "dodgerblue"
            elif nr_pos_eigen == "2":
                col = "orange"
            else:
                col = "crimson"
            ss_steady_state.update({"color": col})
            ss_steady_states += [ss_steady_state]
        return ss_steady_states

    def copy_steady_states(self, b):

        if b.case == 'Trajectories':
            self.customized_points_copy_results.children = []
            results = self.customized_points_copy_results.results
            for index, element in enumerate(self.customized_points_list):
                for key, value in element.iteritems():
                    value.value = results[index][value.description][-1]
        else:
            self.time_courses_copy_final_box.children = []
            self.event_integration_copy_final_box.children = []
            results = b.results
            # Time courses
            for key, value in self.initial_concentration.iteritems():
                value.value = results[value.description][-1]
            # Titrations
            for key, value in self.initial_concentration2.iteritems():
                value.value = results[value.description][-1]

    def export_numerical_results(self, b):

        # Fields of the variable b containing important information are:
        # b.t = contains the time vector
        # b.results = contains a dictionary with relevant data for each variable
        # b.case = contains the type of plot calling the function. Options are: 'Time_Course'; ''; ''

        if b.case == 'Trajectories':

            # create a list of DataFrames
            count = 1
            with pd.ExcelWriter(self.controller.name + '_' + b.case +'.xlsx') as writer:

                for indx, data in enumerate(b.results):
                    order = ['Time'] + data.keys()
                    if self.controller.ds.number_conservations != 0:
                        for n in range(self.controller.ds.number_conservations):
                            if 'Xc' + str(n + 1) in order:
                                order.remove('Xc' + str(n + 1))
                    data.update({'Time': b.t})
                    if self.include_dynamical_in_excel.value is True:
                        data.update({'Dominant Signature': self.dynamical_phenotypes_box.dom_sig[indx],
                                     'Dynamic Phenotype': self.dynamical_phenotypes_box.dynamic_case[indx]})
                        order += ['Dominant Signature', 'Dynamic Phenotype']
                    df = pd.DataFrame(data=data)
                    df = df[order]
                    df.to_excel(writer, sheet_name='Trajectory' + str(count), index=False)
                    count = count + 1

        elif b.case == 'Titrations':

            with pd.ExcelWriter(self.controller.name + '_' + b.case + '.xlsx') as writer:

                order = [b.independent_parameter_identity] + b.results_forward.keys()
                data_forwards = b.results_forward
                df_forwards = pd.DataFrame(data=data_forwards)
                df_forwards.insert(0, b.independent_parameter_identity, b.independent_parameter_forward)
                df_forwards = df_forwards[order]
                df_forwards.to_excel(writer, sheet_name='Forwards', index=False)

                if self.backwards.value is True:
                    data_backwards = b.results_backwards
                    df_backwards = pd.DataFrame(data=data_backwards)
                    df_backwards.insert(0, b.independent_parameter_identity, b.independent_parameter_backwards)
                    df_backwards = df_backwards[order]
                    df_backwards.to_excel(writer, sheet_name='Backwards', index=False)

                el = self.phase_shift_box2
                if el.calculate.value is True and el.target.value != 'None' and el.reference.value != 'None':

                    order = [b.independent_parameter_identity] + [el.target.value]
                    data = b.phase_shift
                    df_response_time = pd.DataFrame(data=data)
                    df_response_time.insert(0, b.independent_parameter_identity, b.independent_parameter_forward)
                    df_response_time = df_response_time[order]
                    df_response_time.to_excel(writer, sheet_name='Phase Shift', index=False)

                    data = b.amplitude
                    df_response_time = pd.DataFrame(data=data)
                    df_response_time.insert(0, b.independent_parameter_identity, b.independent_parameter_forward)
                    df_response_time = df_response_time[order]
                    df_response_time.to_excel(writer, sheet_name='Amplitude', index=False)

                if self.response_times_box.calculate_response_time_titration.value is True:
                    order = [b.independent_parameter_identity] + b.response_time.keys()
                    data = b.response_time
                    df_response_time = pd.DataFrame(data=data)
                    df_response_time.insert(0, b.independent_parameter_identity, b.independent_parameter_forward)
                    df_response_time = df_response_time[order]
                    df_response_time.to_excel(writer, sheet_name='Response Times', index=False)
        else:
            data = b.results
            order = ['Time'] + data.keys()

            if self.controller.ds.number_conservations != 0:
                for n in range(self.controller.ds.number_conservations):
                    if 'Xc' + str(n + 1) in order:
                        order.remove('Xc' + str(n + 1))
            data.update({'Time': b.t})
            if self.include_dynamical_in_excel.value is True:
                data.update({'Dominant Signature': self.dynamical_phenotypes_box.dom_sig,
                             'Dynamic Phenotype': self.dynamical_phenotypes_box.dynamic_case})
                order += ['Dominant Signature', 'Dynamic Phenotype']

            with pd.ExcelWriter(self.controller.name + '_' + b.case +'.xlsx') as writer:

                df1 = pd.DataFrame(data=data)
                df1 = df1[order]
                df1.to_excel(writer, sheet_name='Time Courses', index=False)

                if self.advanced_time_courses.response_time.value != '':
                    response_times = self.time_courses_copy_final_box.response_times
                    order = data.keys()
                    order.remove('Time')
                    df2 = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in response_times.items()]))
                    df2 = df2[order]
                    df2.to_excel(writer, sheet_name='Response Times', index=False)


    def construct_fun(self, y, t):

        # initialize the output array
        y_out = np.array([])

        # Parameter values of this object are actualized to reflect parameter values of the main toolbox.
        for counter, element in enumerate(self.ode_equations):
            self.parameters[str(element.lhs).replace('.', '')] = y[counter]

        # if conserved relationships are present, update the dictionary self.parameters with the respective values.
        if self.conservations_dictionary is not None:
            for conserved_var, expression in self.conservations_dictionary.items():
                self.parameters[conserved_var] = expression.rhs.eval_with_values(p_vals=self.parameters) / self.coef_dic[conserved_var]

        for element in self.ode_equations:
            y_out = np.append(y_out, element.rhs.eval_with_values(p_vals=self.parameters))

        return y_out

    def construct_fun_DAE(self, t, y0, yd0):

        y_out = np.array([])

        # update dictionary of parameter values to consider initial conditions for Xd_t pools
        for indx, element in enumerate(self.ode_equations):
            self.parameters[str(element.lhs).replace('.', '')] = y0[indx]
            self.parameters[str(element.lhs)] = yd0[indx]

        # update dictionary of parameter values to consider intial conditions for Xd_a_c pools
        for indx2, element in enumerate(self.conservations_variable_list):
            self.parameters[element] = y0[indx + indx2 + 1]

        # start building the y_out vector
        for element in self.ode_equations:
            y_out = np.append(y_out, element.rhs.eval_with_values(p_vals=self.parameters))

        for variable in self.conservations_variable_list:
            eq = self.conservations_equal_zero_dictionary[variable]
            y_out = np.append(y_out,
                              eq.rhs.eval_with_values(p_vals=self.parameters))

        y_out = y_out - yd0

        return y_out

    def process_results_vector(self, vector):

        logarithmic = 'Logarithmic' == self.coordinates.value
        if logarithmic is True:
            try:
                processed_vectors = np.log10(vector)
            except:
                print("Divide by zero encountered in log10. Changing to Cartesian coordinates")
                processed_vectors = vector
                self.coordinates.value = "Cartesian"
        else:
            processed_vectors = vector

        return processed_vectors

    def process_label_string(self, s, function=False):
        logarithmic = 'Logarithmic' == self.coordinates.value

        if function is False:
            if s in self.controller.symbols:
                s = '$' + self.controller.symbols[s] + '$'
            processed_s = s if logarithmic is False else r'$\log_{10}$(' + s + ")"
        else:
            processed_s = s if logarithmic is False else 'log(' + s + ")"

        return processed_s

    def make_plot(self, time, results_dic, b):

        logarithmic = 'Logarithmic' == self.coordinates.value
        controller = self.controller
        y = self.process_results_vector(results_dic[self.ylabel.value])
        flag_secondary = False

        # plotting routine
        fig = plt.figure(figsize=[5.5698, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.1714, 0.2, 0.6, 0.7])
        ax.set_title('Time Course Plot')
        ax.plot(time, y, '-')
        ax.set_ylabel(self.process_label_string(self.ylabel.value))

        if self.ylabel_secondary.value != 'None':
            flag_secondary = True
            ax2 = ax.twinx()
            y2 = results_dic[self.ylabel_secondary.value]
            y2 = self.process_results_vector(y2)
            ax2.plot(time, y2, '-', color='orange')
            ax2.set_ylabel(self.process_label_string(self.ylabel_secondary.value))

        for member in self.advanced_time_courses.ss_functions.members:
            if member.name_.value == '' or member.f_.value == '':
                continue
            if member.axis_.value == 'Secondary' and flag_secondary is False:
                flag_secondary = True
                ax2 = ax.twinx()
            if member.axis_.value == "Primary":
                y = self.process_results_vector(results_dic[str(member.name_.value)])
                ax.plot(time, y, color=str(member.col_.value))
                ax.set_ylabel(ax.get_ylabel() + ', ' + self.process_label_string(str(member.name_.value)))
            else:
                y = self.process_results_vector(results_dic[str(member.name_.value)])
                ax2.plot(time, y, color=str(member.col_.value))
                if ax2.get_ylabel() == '':
                    ax2.set_ylabel(self.process_label_string(str(member.name_.value)))
                else:
                    ax2.set_ylabel(ax2.get_ylabel() + ', ' + self.process_label_string(str(member.name_.value)))

        if self.automatic_scale is False:
            ax.set_ylim(self.process_results_vector(self.y_min.value), self.process_results_vector(self.y_max.value))
            if self.ylabel_secondary.value != 'None' or len(self.advanced_time_courses.ss_functions.members) > 0:
                ax2.set_ylim(self.process_results_vector(self.y_min_secondary.value),
                             self.process_results_vector(self.y_max_secondary.value))

        new_caption = ''
        if self.advanced_time_courses.response_time.value != '':

            variable_str = str(self.advanced_time_courses.response_time.value)
            method = self.advanced_time_courses.response_time_method.value
            if method == 'custom' or method == 'y_final*(1-Thr.) ; y_final*(1+Thr.)' or method == 'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ':
                threshold = self.advanced_time_courses.ignore_if_less.value/100
            else:
                threshold = 0
            value = calculate_response_time(time, results_dic, variable_str, method, threshold)
            if len(value) != 0:
                new_caption = self.caption.value + ' The response time for ' +\
                              str(self.advanced_time_courses.response_time.value) + ' is ' + str(value[0]) + '.'

            # Now calculate response time for all other variables.
            response_times_dic = {}
            for variable_str in results_dic.keys():
                # try:
                value = calculate_response_time(time, results_dic, str(variable_str), method, threshold)
                response_times_dic.update({str(variable_str): value})
                # except:
                #     continue
            self.time_courses_copy_final_box.response_times = response_times_dic

        if self.reference_time_course.value != 'None' and self.target_phase_shift_course.value != 'None':
            phase_shift, _ = calculate_phase_shifts(time,
                                                    results_dic,
                                                    self.reference_time_course.value,
                                                    self.target_phase_shift_course.value)
            if len(phase_shift) != 0:
                phase_shift = [str(el) for el in phase_shift]
                if new_caption == '':
                    new_caption = self.caption.value + ' The phase shift values are: ' + ', '.join(phase_shift) + '.'
                else:
                    new_caption += ' The phase shift values are: ' + ', '.join(phase_shift)

        ax.set_xlabel(self.xlabel.value)
        buf = cStringIO.StringIO()
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(buf)
        data = buf.getvalue()
        plt.close()
        controller.figures.add_figure(data,
                                      title=self.title.value,
                                      caption=self.caption.value if new_caption == '' else new_caption,
                                      pvals=b.pvals,
                                      colors=None)

    def make_plot_dynamic_phenotypes(self, time, results_dic, b):

        dynamic_phenotypes = self.dynamical_phenotypes_box.dynamic_case
        controller = self.controller
        y = self.process_results_vector(results_dic[self.ylabel.value])

        # 1. Extract list of unique dynamic phenotypes.
        unique, inverse_indx = np.unique(dynamic_phenotypes, return_inverse=True)
        if len(unique) == 0:
            return

        # 2. Create a dictionary with random colors for each dynamic phenotype
        col_dir = {}
        for el in unique:
            col = "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
            col_dir.update({el: col})

        # 2.5 Initialize plots
        fig = plt.figure(figsize=[5.5698, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.1714, 0.2, 0.6, 0.7])
        ax.set_title('Time Course Plot')

        # 2.6 initialize secondary axis
        if self.ylabel_secondary.value != 'None':
            ax2 = ax.twinx()
            y2 = results_dic[self.ylabel_secondary.value]
            y2 = self.process_results_vector(y2)

        # 3. Loop over dynamic_phenotypes vector and decide when to plot
        masking_vector = np.zeros(len(dynamic_phenotypes), dtype=bool)
        already_plotted = []
        for indx, current_phen in enumerate(dynamic_phenotypes):

            if indx + 1 < len(dynamic_phenotypes): # we are not analyzing the last member of the list

                next_phen = dynamic_phenotypes[indx+1]

                if next_phen != current_phen: # plot and reset masking vector
                    masking_vector[indx] = True
                    time_ = time[masking_vector]
                    y_ = y[masking_vector]

                    if current_phen in already_plotted:
                        ax.plot(time_, y_, '-', color=col_dir[current_phen])
                        ax.set_ylabel(self.process_label_string(self.ylabel.value))
                        if self.ylabel_secondary.value != 'None':
                            y2_ = y2[masking_vector]
                            ax2.plot(time_, y2_, '-', color=col_dir[current_phen])
                            ax2.set_ylabel(self.process_label_string(self.ylabel_secondary.value))
                    else:
                        ax.plot(time_, y_, '-', label='Phen. ' + str(current_phen), color=col_dir[current_phen])
                        ax.set_ylabel(self.process_label_string(self.ylabel.value))
                        if self.ylabel_secondary.value != 'None':
                            y2_ = y2[masking_vector]
                            ax2.plot(time_, y2_, '-', color=col_dir[current_phen])
                            ax2.set_ylabel(self.process_label_string(self.ylabel_secondary.value))
                        already_plotted.append(current_phen)

                    masking_vector = np.zeros(len(dynamic_phenotypes), dtype=bool)

                else: # if the next phenotype is the same, do not plot but update the masking vector
                    masking_vector[indx] = True

            else:

                time_ = time[masking_vector]
                y_ = y[masking_vector]

                if current_phen in already_plotted:
                    ax.plot(time_, y_, '-', color=col_dir[current_phen])
                    ax.set_ylabel(self.process_label_string(self.ylabel.value))
                    if self.ylabel_secondary.value != 'None':
                        y2_ = y2[masking_vector]
                        ax2.plot(time_, y2_, '-', color=col_dir[current_phen])
                        ax2.set_ylabel(self.process_label_string(self.ylabel_secondary.value))
                else:
                    ax.plot(time_, y_, '-', label='Phen. ' + str(current_phen), color=col_dir[current_phen])
                    ax.set_ylabel(self.process_label_string(self.ylabel.value))
                    if self.ylabel_secondary.value != 'None':
                        y2_ = y2[masking_vector]
                        ax2.plot(time_, y2_, '-', color=col_dir[current_phen])
                        ax2.set_ylabel(self.process_label_string(self.ylabel_secondary.value))
                    already_plotted.append(current_phen)

                masking_vector = np.zeros(len(dynamic_phenotypes), dtype=bool)

        if self.automatic_scale is False:
            ax.set_ylim(self.process_results_vector(self.y_min.value), self.process_results_vector(self.y_max.value))
            if self.ylabel_secondary.value != 'None' or len(self.advanced_time_courses.ss_functions.members) > 0:
                ax2.set_ylim(self.process_results_vector(self.y_min_secondary.value),
                             self.process_results_vector(self.y_max_secondary.value))

        ax.set_xlabel(self.xlabel.value)
        ax.legend(fontsize=5, loc='best', frameon=False)
        buf = cStringIO.StringIO()
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(buf)
        data = buf.getvalue()
        plt.close()
        controller.figures.add_figure(data,
                                      title=self.title.value,
                                      caption=self.caption.value,
                                      pvals=b.pvals,
                                      colors=None)

    def draw_1D_ss_function(self, ax, function, p_vals,
                            slice_variable, range_slice,
                            resolution=100, colors=None, included_cases=None, ylim=None, **kwargs):
        object = self.controller.ds
        p_bounds = dict(p_vals)
        p_bounds[slice_variable] = range_slice
        if included_cases is not None:
            included_cases = [i.case_number for i in object(included_cases)]
            if object.number_of_cases < 1e5:
                valid_cases = object.valid_cases(p_bounds=p_bounds)
                valid_nonstrict = boject.valid_cases(p_bounds=p_bounds, strict=False)
                hatched_cases = [i for i in valid_cases if i not in included_cases]
                valid_nonstrict = [i for i in valid_cases if i in included_cases]
                valid_cases = [i for i in valid_cases if i in included_cases]
                valid_nonstrict = [i for i in valid_cases if i not in valid_cases]
            else:
                valid_cases = [i for i in included_cases if object(i).is_valid(p_bounds=p_bounds)]
                valid_nonstrict = []
        else:
            valid_cases = object.valid_cases(p_bounds=p_bounds)
            valid_nonstrict = object.valid_cases(p_bounds=p_bounds, strict=False)
            valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
        if len(valid_cases) + len(valid_nonstrict) == 0:
            # fill black
            return
        valid_cases = valid_cases + valid_nonstrict  # self.valid_cases(p_bounds=p_bounds)
        lines = list()
        ylim_t = None
        if 'color' not in kwargs:
            kwargs['color'] = 'k'
        for case in valid_cases:
            if colors is not None:
                if case in colors:
                    kwargs['color'] = colors[case]
                else:
                    kwargs.pop('color')
            pt = object(case).draw_1D_ss_function(ax,
                                                function,
                                                p_vals,
                                                slice_variable,
                                                range_slice,
                                                resolution=resolution,
                                                **kwargs)
            if pt is None:
                continue
            lines.append(pt)
            ydata = pt[0].get_ydata()
            miny = min(ydata)
            maxy = max(ydata)
            if ylim is None:
                if ylim_t is None:
                    ylim_t = [min(ydata), max(ydata)]
                else:
                    ylim_t = [min((ylim_t[0], miny)), max((ylim_t[1], maxy))]
        if ylim is None:
            ylim = ylim_t
        ax.set_ylim([ylim[0] - (ylim[1] - ylim[0]) * 0.1, ylim[1] + (ylim[1] - ylim[0]) * 0.1])
        ax.set_xlim(np.log10(range_slice))
        if slice_variable in object._latex:
            slice_variable = '$' + object._latex[slice_variable] + '$'
        ax.set_xlabel(r'$\log_{10}$(' + slice_variable + ')')
        if isinstance(function, Expression):
            expr = function
        else:
            expr = Expression(function)
        # ax.set_ylabel('$' + expr.latex(object._latex) + '$')
        return lines


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)


def find_crossings(y, y_target):
    y_bol = [y_target < el for el in y ]
    y_bol
    crossings = []
    search_crossing = True
    for i, element in enumerate(y_bol):
        if search_crossing is True:
            bool_value = element
            search_crossing = False
            continue
        if element != bool_value:
            if bool_value is True and element is False:
                crossings.append([i-1, i, 'Decreasing'])
                bool_value = element
            else:
                crossings.append([i-1, i, 'Increasing'])
                bool_value = element
    return crossings


def interpolate_crossings(t, y, crossings, y_target):
    value = []
    for crossing in crossings:
        x1 = t[crossing[0]]
        y1 = y[crossing[0]]
        x2 = t[crossing[1]]
        y2 = y[crossing[1]]
        m = (y1 - y2)/(x1 - x2)
        x_target = (y_target - y1)/m + x1
        value.append(round_sig(x_target))
    return value


def calculate_response_time(time, results_dic, variable_str, method, threshold):

    ## The different methods include: abs(y_final - y_initial); max(y) - min(y); custom
    y = results_dic[str(variable_str)]
    value = calculate_response_time_vector(time, y, method, threshold)

    return value

def calculate_response_time_vector(time, y, method, threshold):
# def calculate_response_time_vector(time, results_dic, variable_str, method, threshold):
    ## The different methods include: abs(y_final - y_initial); max(y) - min(y); custom

    value = []

    if method == 'y_final - y_initial':
        # y = results_dic[str(variable_str)]
        d = abs(abs(y[0]) - abs(y[-1])) / 2
        y_target = d + y[0] if y[0] < y[-1] else d + y[-1]
        crossings = find_crossings(y, y_target)
        value = interpolate_crossings(time, y, crossings, y_target)

    elif method == 'max(y) - min(y)':
        # y = results_dic[str(variable_str)]
        d = abs(abs(max(y)) - abs(min(y))) / 2
        y_target = max(y) - d
        crossings = find_crossings(y, y_target)
        value = interpolate_crossings(time, y, crossings, y_target)

    elif method == 'custom':

        # y = results_dic[str(variable_str)]

        # The first case describes the situation where the initial point is the minimum.
        #  The time course goes up and then down.

        if (y[0] == min(y) or y[0] > y[-1]) and max(y) != y[0]:

            # catch the first response time. Going up
            d = abs(abs(y[0]) - abs(max(y))) / 2
            y_target = y[0] + d
            crossings = find_crossings(y, y_target)
            value1 = interpolate_crossings(time, y, crossings, y_target)
            if len(value1) > 0:
                value.append(value1[0])

            # catch the second response time, when it goes down
            if max(y) - y[-1] > threshold * (max(y) - min(y)):
                d = abs(abs(max(y)) - abs(y[-1])) / 2
                y_target = max(y) - d
                crossings = find_crossings(y, y_target)
                value2 = interpolate_crossings(time, y, crossings, y_target)
                if len(value2) > 1:
                    value.append(value2[1])
                else:
                    if len(value2) > 0:
                        value.append(value2[0])

        # The second case describes the situation where the initial point is the maximum.
        #  The time course goes down and then up.

        else:

            # finding the response time when curve goes down.
            if y[0] > min(y):
                d = abs(abs(y[0]) - abs(min(y))) / 2
                y_target = y[0] - d
                crossings = find_crossings(y, y_target)
                value1 = interpolate_crossings(time, y, crossings, y_target)
                if len(value1) > 0:
                    value.append(value1[0])

            # Catch the second response time when the curve goes up
            if y[-1] - min(y) > threshold*(max(y) - min(y)):
                d = abs(abs(y[-1]) - abs(min(y))) / 2
                y_target = min(y) + d
                crossings = find_crossings(y, y_target)
                value2 = interpolate_crossings(time, y, crossings, y_target)
                if len(value2) > 1:
                    value.append(value2[1])
                else:
                    if len(value2) > 0:
                        value.append(value2[0])

        if len(value) == 1:
            value.append(value[0])

    elif method == 'y_final*(1-Thr.) ; y_final*(1+Thr.)':

        y_up = y[-1] * (1 + threshold)
        y_down = y[-1] * (1 - threshold)

        crossings_up = find_crossings(y, y_up)
        value_up = interpolate_crossings(time, y, crossings_up, y_up)
        if len(value_up) != 0:
            value_up = value_up[-1]
        else:
            value_up = 0

        crossings_down = find_crossings(y, y_down)
        value_down = interpolate_crossings(time, y, crossings_down, y_down)
        if len(value_down) != 0:
            value_down = value_down[-1]
        else:
            value_down = 0

        value = [value_up] if value_up > value_down else [value_down]

    elif method == 'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ':

        y_up = y[-1] + abs(y[-1] - y[0]) * threshold
        y_down = y[-1] - abs(y[-1] - y[0]) * threshold

        crossings_up = find_crossings(y, y_up)
        value_up = interpolate_crossings(time, y, crossings_up, y_up)
        if len(value_up) != 0:
            value_up = value_up[-1]
        else:
            value_up = 0

        crossings_down = find_crossings(y, y_down)
        value_down = interpolate_crossings(time, y, crossings_down, y_down)
        if len(value_down) != 0:
            value_down = value_down[-1]
        else:
            value_down = 0

        value = [value_up] if value_up > value_down else [value_down]

    return value


def update_if_less_function(name, value, widget=None):
    if value == 'custom' or value == 'y_final*(1-Thr.) ; y_final*(1+Thr.)' or value == 'y_final + abs(y_final - y_initial)*Thr. ; y_final - abs(y_final - y_initial)*Thr. ' and widget is not None:
        widget.visible = True
    else:
        widget.visible = False


def calculate_phase_shifts(time, results_dic, reference_time_course, target_phase_shift_course, mode='legend'):

    # Get time courses for reference and target values.
    reference_vector = results_dic[reference_time_course]
    target_vector = results_dic[target_phase_shift_course]
    phase_shift, amplitude = calculate_phase_shifts_from_vectors(time, reference_vector, target_vector, mode=mode)

    return phase_shift, amplitude


def calculate_phase_shifts_from_vectors(time, reference_vector, target_vector, mode='legend'):

    # Get peaks
    ref_peaks = detect_peaks(reference_vector)
    target_peaks = detect_peaks(target_vector)

    # Get valleys for target
    target_valleys = detect_peaks(target_vector, valley=True)

    if mode == 'legend':

        # For each peak in ref_vals, calculate phase shift.
        phase_shift = set()
        for peak in ref_peaks:
            previous_target_peaks = target_peaks[target_peaks < peak]
            if len(previous_target_peaks) != 0:
                target_peak = previous_target_peaks[-1]
                value = round_sig(time[peak] - time[target_peak])
                phase_shift.add(value)

        # For each peak in target_peaks, calculate amplitude
        amplitude = set()
        for peak in target_peaks:
            previous_target_valleys = target_valleys[target_valleys < peak]
            if len(previous_target_valleys) != 0:
                valley = previous_target_valleys[-1]
                value = round_sig(target_vector[peak] - target_vector[valley])
                amplitude.add(value)
    else:

        # For each peak in ref_vals, calculate phase shift.
        phase_shift = np.array([])
        for peak in ref_peaks:
            previous_target_peaks = target_peaks[target_peaks < peak]
            if len(previous_target_peaks) != 0:
                target_peak = previous_target_peaks[-1]
                value = round_sig(time[peak] - time[target_peak])
                phase_shift = np.append(phase_shift, value)
        phase_shift = phase_shift.mean() if len(phase_shift) > 0 else np.nan

        # For each peak in target_peaks, calculate amplitude
        amplitude = np.array([])
        for peak in target_peaks:
            previous_target_valleys = target_valleys[target_valleys < peak]
            if len(previous_target_valleys) != 0:
                valley = previous_target_valleys[-1]
                value = round_sig(target_vector[peak] - target_vector[valley])
                amplitude = np.append(amplitude, value)
        amplitude = amplitude.mean() if len(amplitude) > 0 else np.nan



    return phase_shift, amplitude


def find_vector(variable_list, y, variable):
    vector = []
    for counter, element in enumerate(variable_list):
        if element == variable:
            vector = y[:, counter]
    return vector

