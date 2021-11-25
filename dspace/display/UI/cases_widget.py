from __future__ import division
import dspace
import dspace.plotutils
import dspace.display
import time
from dspace.SWIG.dspace_interface import DSCaseFree
from decimal import Decimal
from math import log10, floor
import pandas as pd
from dspace.SWIG.dspace_interface import DSDesignSpaceNumberOfBoundaries


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



    Select = RadioButtons
    VBox = Box
    HBox = Box
    ipy_old = True
else:
    from ipywidgets import *
    from popup import Popup

    ipy_old = False

from IPython.display import clear_output, display, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg

from case_widget import DisplayCase


class CasesTable(object):

    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'table', VBox())

    def cases_table_widget(self):
        controller = self.controller
        options = Select(values=[#'None',
                                 'All',
                                 'Valid',
                                 'User Specified'],
                         options=[#'None',
                                  'All',
                                  'Valid',
                                  'User Specified'],
                         value='Valid')
        if 'biological_constraints' not in controller.options:
            bio_constraints = ''
        else:
            bio_constraints = ', '.join(controller.defaults('biological_constraints'))
        cases = Textarea(description='User specified cases')
        by_signature = Checkbox(description='Cases indicated by signature?', value=False)
        constraints = Textarea(description='Biological constraints:',
                               value=bio_constraints)
        show_signature = Checkbox(description='Show signature?', value=True)
        show_parents = Checkbox(description='Show parent cases?', value=False)
        export_table = Checkbox(description='Export Phenotypic Repertoire (.xlsx)', value=False)
        add_column = Button(description='Add column')
        add_filter = Button(description='Add column filter')
        remove_column = Button(description='Remove column')
        remove_column.visible = False
        add_column.remove_column = remove_column
        extra_columns = VBox(children=[])
        self.extra_columns = extra_columns
        filters = VBox(children=[])
        add_column.on_click(self.add_cases_column)
        add_filter.on_click(self.add_cases_filter)
        add_column.column_block = extra_columns
        add_filter.column_block = extra_columns
        add_filter.filters = filters
        remove_column.on_click(self.remove_cases_column)
        remove_column.column_block = extra_columns
        case_id = Text(description='Report for case:')
        wi = VBox(children=[options,
                            cases,
                            by_signature,
                            constraints,
                            show_signature, show_parents, export_table,
                            add_column,
                            extra_columns,
                            remove_column,
                            filters,
                            add_filter
                            ])
        button = Button(value=False, description='Enumerate Phenotypic Repertoire')  #'Create/modify cases table'
        button.on_click(self.show_cases)
        button.wi = wi
        button.options = options
        button.cases = cases
        button.by_signature = by_signature
        button.extra_columns = extra_columns
        button.constraints = constraints
        button.filters = filters
        button.show_signature = show_signature
        button.show_parents = show_parents
        button.export_table = export_table
        cases_table = VBox(description='Cases Table', children=[wi, button, self.table])
        wi.visible = False
        return ('Phenotypic Repertoire', cases_table)

    def add_cases_filter(self, b):
        # if len(b.column_block.children) == 0:
        #     return
        controller = self.controller
        children = [i for i in b.filters.children]
        new_filter = [Dropdown(description='Column #',
                               values=[str(i + 1) for i in range(1, 2+len(b.column_block.children))],
                               options=[str(i + 1) for i in range(1, 2+len(b.column_block.children))],
                               value=str(2+len(b.column_block.children))),
                      Dropdown(description='Condition',
                               values=['==', '>', '<', '>=', '<=', '!='],
                               options=['==', '>', '<', '>=', '<=', '!=']),
                      Text(description='Value'),
                      Button(description='X')]
        new_filter[3].on_click(self.remove_cases_filter)
        current = HBox(children=new_filter)
        new_filter[3].current = current
        new_filter[3].filters = b.filters
        children.append(current)
        b.filters.children = children

    # def update_cases_filter(self, filters, columns):

    def remove_cases_filter(self, b):
        controller = self.controller
        if len(b.filters.children) == 0:
            return
        children = [i for i in b.filters.children]
        for i, child in enumerate(children):
            if b.current == child:
                children.pop(i)
        b.filters.children = children

    def update_column_block(self, name, value):

        for column in self.extra_columns.children:

            if column.header.value == 'Validity' or column.header.value == 'Is a Blowup':
                column.dependent.visible = False
                column.independent.visible = False
                # column.ignore_infinity.visible = False
                # column.log_coordinates.visible = False
                column.root_method.visible = False
                column.show_complex_conjugates.visible = False
                column.volume_method.visible = False
                column.nr_vertices.visible = False
                column.lower_bound_par_space.visible = False
                column.upper_bound_par_space.visible = False
                column.limit_vertices.visible = False

            elif column.header.value == 'Volume' or column.header.value == 'Normalized Volume':
                column.dependent.visible = False
                column.independent.visible = False
                column.root_method.visible = False
                column.show_complex_conjugates.visible = False
                column.volume_method.visible = True
                column.lower_bound_par_space.visible = True
                column.upper_bound_par_space.visible = True
                column.limit_vertices.visible = False

                if column.volume_method.value == 'Tolerances' or \
                        column.volume_method.value == 'Bounding Box' or column.volume_method.value == 'Geometric Mean T. & BB.'  :
                    # column.ignore_infinity.visible = True
                    # column.log_coordinates.visible = True
                    column.nr_vertices.visible = False
                    column.limit_vertices.visible = False
                else:  # lrs library
                    # column.ignore_infinity.visible = False
                    # column.log_coordinates.visible = False
                    column.lower_bound_par_space.visible = True
                    column.upper_bound_par_space.visible = True
                    column.limit_vertices.visible = True
                    if column.limit_vertices.value is True:
                        column.nr_vertices.visible = True
                    else:
                        column.nr_vertices.visible = False
            elif column.header.value == 'Log-Gain':
                column.dependent.visible = True
                column.independent.visible = True
                # column.ignore_infinity.visible = False
                # column.log_coordinates.visible = False
                column.root_method.visible = False
                column.show_complex_conjugates.visible = False
                column.volume_method.visible = False
                column.lower_bound_par_space.visible = False
                column.upper_bound_par_space.visible = False
                column.limit_vertices.visible = False

            elif column.header.value == 'Dimension':
                column.dependent.visible = False
                column.independent.visible = False
                column.root_method.visible = False
                column.show_complex_conjugates.visible = False
                column.volume_method.visible = False
                column.nr_vertices.visible = False
                # column.ignore_infinity.visible = False
                # column.log_coordinates.visible = False
                column.lower_bound_par_space.visible = True
                column.upper_bound_par_space.visible = True
                column.limit_vertices.visible = False

            else:    # column.header.value == '# eigenvalues w/ positive real part':
                column.dependent.visible = False
                column.independent.visible = False
                # column.ignore_infinity.visible = False
                # column.log_coordinates.visible = False
                column.root_method.visible = True
                column.volume_method.visible = False
                column.limit_vertices.visible = False

                if column.root_method.value == 'Numpy Functions':
                    column.show_complex_conjugates.visible = True
                if column.root_method.value == 'Routh Array':
                    column.show_complex_conjugates.visible = False

        column.ignore_infinity.visible = False
        column.log_coordinates.visible = False

    def add_cases_column(self, b):
        controller = self.controller
        nr_ineq = DSDesignSpaceNumberOfBoundaries(controller.ds._swigwrapper) + \
                  len(controller.constraints) + len(controller.pvals)*2
        try:
            max_vertices = nr_ineq ** (len(controller.pvals) / 2)
        except:
            max_vertices = 1e300

        children = [i for i in b.column_block.children]
        new_column = [Dropdown(description='Column ' + str(len(children) + 3),                          #0
                               values=['Validity',
                                       'Is a Blowup',
                                       'Volume',
                                       'Normalized Volume',
                                       '# eigenvalues w/ positive real part',
                                       'Log-Gain',
                                       'Dimension'],
                               options=['Validity',
                                        'Is a Blowup',
                                        'Volume',
                                        'Normalized Volume',
                                        '# eigenvalues w/ positive real part',
                                        'Log-Gain',
                                        'Dimension'],
                               value='Log-Gain'),
                      Dropdown(description='Log-gain (dependent variable)',                             #1
                               values=controller.ds.dependent_variables,
                               options=controller.ds.dependent_variables),

                      Dropdown(description='Log-gain (independent variable)',                           #2
                               values=controller.ds.independent_variables,
                               options=controller.ds.independent_variables),

                      Dropdown(description='Volume Calculation Method',                                 #3
                               values=['Vertex Enumeration',
                                       'Tolerances',
                                       'Bounding Box',
                                       'Geometric Mean T. & BB.',
                                       ],
                               options=['Vertex Enumeration',
                                        'Tolerances',
                                        'Bounding Box',
                                        'Geometric Mean T. & BB.'
                                        ],
                               value='Tolerances',
                               visible=False),

                      Checkbox(description='Limit Number of Vertices',                                  #4
                               value=False,
                               visible=False),

                      Slider(description='Max. Number of Vertices',                                      #5
                             value=int(max_vertices),
                             min=1,
                             max=round_sig(max_vertices, 1)*2,
                             step=round_sig(max_vertices, 1)/100,
                             visible=False),

                      Slider(description='Upper Boundary Par. Space (Log)',                              #6
                             value=3,
                             min=1,
                             max=20,
                             step=0.25,
                             visible=False),

                      Slider(description='Lower Boundary Par. Space (Log)',                              #7
                             min=-20,
                             max=-1,
                             step=0.25,
                             value=-3,
                             visible=False),

                      Dropdown(description='Ignore Unbounded Parameters',                               #8
                               values=[True, False],
                               options=['True', 'False'],
                               value=False,
                               visible=False),
                      Dropdown(description='Logarithmic Coordinates',                                   #9
                               values=[True, False],
                               options=['True', 'False'],
                               visible=False),
                      Dropdown(description='Calculate Positive Roots Using:',                           #10
                               values=['Routh Array', 'Numpy Functions'],
                               options=['Routh Array', 'Numpy Functions'],
                               value='Routh Array',
                               visible=False),
                      Dropdown(description='Check for complex conjugates:',                             #11
                               values=[True, False],
                               options=['True', 'False'],
                               value=True,
                               visible=False),
                      ]

        column = VBox(children=new_column)
        column.header = new_column[0]
        column.dependent = new_column[1]
        column.independent = new_column[2]
        column.volume_method = new_column[3]

        column.limit_vertices = new_column[4]

        column.nr_vertices = new_column[5]

        column.upper_bound_par_space = new_column[6]
        column.lower_bound_par_space = new_column[7]

        column.ignore_infinity = new_column[8]
        column.log_coordinates = new_column[9]
        column.root_method = new_column[10]
        column.show_complex_conjugates = new_column[11]
        children.append(column)
        b.column_block.children = children
        column.header.on_trait_change(self.update_column_block, 'value')
        column.root_method.on_trait_change(self.update_column_block, 'value')
        column.volume_method.on_trait_change(self.update_column_block, 'value')
        column.limit_vertices.on_trait_change(self.update_column_block, 'value')
        b.remove_column.visible = True

    def remove_cases_column(self, b):
        controller = self.controller
        if len(b.column_block.children) == 0:
            return
        if len(b.column_block.children) == 1:
            b.visible = False
        children = [i for i in b.column_block.children[:-1]]
        b.column_block.children = children

    def show_cases(self, b):
        controller = self.controller
        if controller.ds == None:
            b.wi.visible = False
            self.table.children = []
            b.description = 'Enumerate Phenotypic Repertoire' #  'Create/modify cases table'
            return
        if b.wi.visible == False:
            self.table.children = []
            b.wi.visible = True
            b.description = 'Done'
            return
        self.create_case_table(b, mode=str(b.options.value))
        b.wi.visible = False
        b.description = 'Enumerate Phenotypic Repertoire' #  'Create/modify cases table'

    def create_string_for_table_long(self, cases, b, s):
        d = dict((i, []) for i in b.headers)
        controller = self.controller
        parents = []
        if controller.ds.instability is True:
            if b.show_signature.value is True:
                for c in cases:
                    # if (c.is_valid() is False or c.is_valid(p_bounds=c.valid_interior_parameter_set()) is False) and b.mode == "Valid":
                    if c.is_valid() is False and b.mode == "Valid":
                        continue
                    if b.show_parents.value is True:
                        p = c.case_number.split('_')
                        if len(p) - 1 != 0:
                            case_number = '_'.join([p[i] for i in range(0, len(p)-1)])
                            if case_number not in parents:
                                parents.append(case_number)
                                parent = self.controller.ds(case_number)
                                values = [parent.case_number, parent.signature]
                                values += [self.value_for_extra_column(parent, column) for column in
                                           b.extra_columns.children]
                                if self.show_case(values, b.extra_columns.children, b.filters) is False:
                                    continue
                                if b.export_table.value is True:
                                    for i in range(len(b.headers)):
                                            d[b.headers[i]].append(values[i])
                                s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>' \
                                     '<td style="padding:0 15px 0 15px;">{1}</td>'.format(
                                    parent.case_number + '&#1007' if (
                                                parent.is_cyclical is False and parent.is_unstable is True) else parent.case_number,
                                    parent.signature)
                                for column in b.extra_columns.children:
                                    s += self.html_for_extra_column(parent, column)
                                s += '</tr>\n'
                    values = [c.case_number, c.signature]
                    values += [self.value_for_extra_column(c, column) for column in b.extra_columns.children]
                    if self.show_case(values, b.extra_columns.children, b.filters) is False:
                        continue
                    if b.export_table.value is True:
                        for i in range(len(b.headers)):
                                d[b.headers[i]].append(values[i])
                    s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td><td style="padding:0 15px 0 15px;">{1}</td>'.format(
                         c.case_number+'&#1007' if (c.is_cyclical is False and c.is_unstable is True) else c.case_number, c.signature)
                    for column in b.extra_columns.children:
                        s += self.html_for_extra_column(c, column)
                    s += '</tr>\n'

            else:
                for c in cases:
                    # if (c.is_valid() is False or c.is_valid(p_bounds=c.valid_interior_parameter_set()) is False) and b.mode == "Valid":
                    if c.is_valid() is False and b.mode == "Valid":
                        continue
                    if b.show_parents.value is True:
                        p = c.case_number.split('_')
                        if len(p) - 1 != 0:
                            case_number = '_'.join([p[i] for i in range(0, len(p)-1)])
                            if case_number not in parents:
                                parents.append(case_number)
                                parent = self.controller.ds(case_number)
                                values = [parent.case_number]
                                values += [self.value_for_extra_column(parent, column) for column in
                                           b.extra_columns.children]
                                if self.show_case(values, b.extra_columns.children, b.filters) is False:
                                    continue
                                if b.export_table.value is True:
                                    for i in range(len(b.headers)):
                                            d[b.headers[i]].append(values[i])
                                s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>'.format(
                                    parent.case_number + '&#1007' if (parent.is_cyclical is False and
                                                                      parent.is_unstable is True) else parent.case_number)
                                for column in b.extra_columns.children:
                                    s += self.html_for_extra_column(parent, column)
                                s += '</tr>\n'
                    values = [c.case_number]  # [c.case_number, c.signature]
                    values += [self.value_for_extra_column(c, column) for column in b.extra_columns.children]
                    if self.show_case(values, b.extra_columns.children, b.filters) is False:
                        continue
                    if b.export_table.value is True:
                        for i in range(len(b.headers)):
                                d[b.headers[i]].append(values[i])
                    s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>'.format(
                         c.case_number+'&#1007' if (c.is_cyclical is False and c.is_unstable is True) else c.case_number)
                    for column in b.extra_columns.children:
                        s += self.html_for_extra_column(c, column)
                    s += '</tr>\n'
        else:
            if b.show_signature.value is True:
                for c in cases:
                    # if (c.is_valid() is False or c.is_valid(p_bounds=c.valid_interior_parameter_set()) is False) and b.mode == "Valid":
                    if c.is_valid() is False and b.mode == "Valid":
                        continue
                    if b.show_parents.value is True:
                        p = c.case_number.split('_')
                        if len(p) - 1 != 0:
                            case_number = '_'.join([p[i] for i in range(0, len(p)-1)])
                            if case_number not in parents:
                                parents.append(case_number)
                                parent = self.controller.ds(case_number)
                                values = [parent.case_number, parent.signature]
                                values += [self.value_for_extra_column(parent, column) for column in
                                           b.extra_columns.children]
                                if self.show_case(values, b.extra_columns.children, b.filters) is False:
                                    continue
                                if b.export_table.value is True:
                                    for i in range(len(b.headers)):
                                            d[b.headers[i]].append(values[i])
                                s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>' \
                                     '<td style="padding:0 15px 0 15px;">{1}</td>'.format(
                                    parent.case_number + '&#1007' if (
                                                parent.is_cyclical is False and parent.is_unstable is True) else parent.case_number,
                                    parent.signature)
                                for column in b.extra_columns.children:
                                    s += self.html_for_extra_column(parent, column)
                                s += '</tr>\n'
                    # self.show_case(c, b.extra_columns)
                    values = [c.case_number, c.signature]
                    values += [self.value_for_extra_column(c, column) for column in b.extra_columns.children]
                    if self.show_case(values, b.extra_columns.children, b.filters) is False:
                        continue
                    if b.export_table.value is True:
                        for i in range(len(b.headers)):
                                d[b.headers[i]].append(values[i])
                    s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td><td style="padding:0 15px 0 15px;">{1}</td>'.format(
                         c.case_number, c.signature)
                    for column in b.extra_columns.children:
                        s += self.html_for_extra_column(c, column)
                    s += '</tr>\n'
            else:
                for c in cases:
                    # if (c.is_valid() is False or c.is_valid(p_bounds=c.valid_interior_parameter_set()) is False) and b.mode == "Valid":
                    if c.is_valid() is False and b.mode == "Valid":
                        continue
                    if b.show_parents.value is True:
                        p = c.case_number.split('_')
                        if len(p) - 1 != 0:
                            case_number = '_'.join([p[i] for i in range(0,len(p)-1)])
                            if case_number not in parents:
                                parents.append(case_number)
                                parent = self.controller.ds(case_number)
                                values = [parent.case_number]
                                values += [self.value_for_extra_column(parent, column) for column in
                                           b.extra_columns.children]
                                if self.show_case(values, b.extra_columns.children, b.filters) is False:
                                    continue
                                if b.export_table.value is True:
                                    for i in range(len(b.headers)):
                                            d[b.headers[i]].append(values[i])
                                s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>'.format(
                                    parent.case_number + '&#1007' if (parent.is_cyclical is False and parent.is_unstable
                                                                      is True) else parent.case_number)
                                for column in b.extra_columns.children:
                                    s += self.html_for_extra_column(parent, column)
                                s += '</tr>\n'
                    # self.show_case(c, b.extra_columns)
                    values = [c.case_number] # [c.case_number, c.signature]
                    values += [self.value_for_extra_column(c, column) for column in b.extra_columns.children]
                    if self.show_case(values, b.extra_columns.children, b.filters) is False:
                        continue
                    if b.export_table.value is True:
                        for i in range(len(b.headers)):
                            d[b.headers[i]].append(values[i])
                    s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>'.format(c.case_number)
                    for column in b.extra_columns.children:
                        s += self.html_for_extra_column(c, column)
                    s += '</tr>\n'
        b.table_data = d
        return s

    def create_string_for_table(self, b, cases, s):

        for c in cases:
            values = [c.case_number, c.signature]
            values += [self.value_for_extra_column(c, column) for column in b.extra_columns.children]
            if self.show_case(values, b.extra_columns.children, b.filters) is False:
                continue
            s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td><td style="padding:0 15px 0 15px;">{1}</td>'.format(
                c.case_number, c.signature)
            for column in b.extra_columns.children:
                s += self.html_for_extra_column(c, column)
            s += '</tr>\n'
        return s

    def create_case_table(self, b, mode='All'):
        controller = self.controller
        constraints = str(b.constraints.value)
        headers = []
        b.mode = mode
        if constraints != '':
            constraints = constraints.split(',')
            constraints = [i.strip() for i in constraints if len(i.strip()) > 0]
        else:
            constraints = []
        controller.set_defaults('biological_constraints', constraints)
        if mode == 'None':
            self.table.children = []
            return
        s = '<div><table>\n<caption>Cases in the system design space. </caption>\n'
        if b.show_signature.value:
            s += '<tr align=center><td style="padding:0 15px 0 15px;"><b>{0}</b></td><td style="padding:0 15px 0 15px;">' \
                 '<b>{1}</b></td>'.format(
                 '  Case Number  ', '  Case Signature  ')
            headers += ['Case Number', 'Case Signature']
        else:
            s += '<tr align=center><td style="padding:0 15px 0 15px;"><b>{0}</b></td>'.format('  Case Number  ')
            headers += ['Case Number']

        for column in b.extra_columns.children:
            header = str(column.header.value)
            header_log = 'Log ' if column.log_coordinates.value is True else ''

            if header == 'Log-Gain':
                header = 'L(' + str(column.dependent.value) + ',' + str(column.independent.value) + ')'

            elif header == 'Is a Blowup':
                header = 'Is a Blowup'

            elif header == 'Volume':
                if column.volume_method.value == 'Tolerances':
                    header = 'Volume <br>' + '(ignoring unbnd prmtrs; Tolerances)' if column.ignore_infinity.value is True \
                        else 'Volume <br>(Tolerances)'
                    header = header_log + header
                elif column.volume_method.value == 'Bounding Box':
                    header = 'Volume <br>' + '(ignoring unbnd prmtrs; Bounding Box)' if column.ignore_infinity.value is True \
                        else 'Volume <br>(Bounding Box)'
                    header = header_log + header
                elif column.volume_method.value == 'Geometric Mean T. & BB.':
                    header = 'Volume <br> (Geometric Mean)'
                    header = header_log + header
                else:
                    header = 'Log Volume<br>(lrs)'

            elif header == 'Dimension':
                header = 'Dimension'

            elif header == 'Normalized Volume':
                method = column.volume_method.value
                if header_log == 'Log ':
                    header = 'Normalized Log Volume <br>' + '(ignoring unbnd prmtrs; ' + method + ')' if \
                        column.ignore_infinity.value is True else 'Normalized Log Volume <br>' + '(' + method + ')'
                else:
                    header = 'Normalized Volume <br>' + '(ignoring unbnd prmtrs; ' + method + ')' if \
                        column.ignore_infinity.value is True else 'Normalized Volume <br>' + '(' + method + ')'

            elif header == '# eigenvalues w/ positive real part':
                if column.root_method.value == 'Numpy Functions':
                    header = '# eigenvalues w/ positive real part (Numpy)'
                    if column.show_complex_conjugates.value is True:
                        header = 'Has complex conjugates'

            s += '<td><b>' + header + '</b></td>'
            headers += [header]
        s += '</tr>'
        if mode == 'Valid':
            cases = controller.ds(controller.ds.valid_cases(),
                                  constraints=constraints)
        elif mode == 'All':
            cases = controller.ds(
                [i for i in range(1, controller.ds.number_of_cases + 1)],
                constraints=constraints)
        else:
            case_signatures = [i.strip() for i in str(b.cases.value).split(',') if len(i.strip()) > 0]
            cases = controller.ds(["".join(i.split()) for i in case_signatures],
                                  by_signature=b.by_signature.value,
                                  constraints=constraints)

        if mode == 'Valid' or mode == 'User Specified':
            self.calculate_total_volumes(cases, b, log_coordinate=False)

        b.headers = headers
        s = self.create_string_for_table_long(cases, b, s)
        if b.export_table.value is True:
            df = pd.DataFrame(data=b.table_data)
            df.to_excel(controller.name + '.xlsx', index=False, columns=headers)
        # s = self.create_string_for_table(cases, b, s)

        s += '</table><caption>'
        s += 'Note: # of eigenvalues w/ positive real part is calculated using a representative set of parameter values, and may not be reflective of all potential behaviors.'
        s += '<br>The number of cases listed in this table is: ' + str(s.count("</tr>\n"))
        s += '</caption></div>'

        html_widget = HTML(value=s)
        save_table = Button(description='Save Table')
        save_table.table_data = s
        save_table.on_click(self.save_table)
        table_container = VBox(children=[save_table,
                                         html_widget,
                                         ])
        if ipy_old is True:
            table_container.set_css('height', '300px')
        else:
            table_container.height = '300px'
            table_container.overflow_x = 'auto'
            table_container.overflow_y = 'auto'
        self.table.children = [table_container]

    def calculate_total_volumes(self, cases, b, log_coordinate=True):

        n = 0
        for column in b.extra_columns.children:

            method = column.volume_method.value
            ignore_unbounded = column.ignore_infinity.value
            log_coordinate = column.log_coordinates.value
            lowerBounds = 10 ** column.lower_bound_par_space.value
            upperBounds = 10 ** column.upper_bound_par_space.value
            maxVertices = column.nr_vertices.value
            limitVertices = column.limit_vertices.value

            if column.header.value == 'Normalized Volume':
                vol = 0
                for c in cases:
                    if c.is_valid() is False or c.is_valid(p_bounds=c.valid_interior_parameter_set()) is False:
                        continue
                    values = [c.case_number, c.signature]
                    values += [self.value_for_extra_column(c, column, filtering_mode=True) for column in b.extra_columns.children]
                    if self.show_case(values, b.extra_columns.children, b.filters) is False:
                        continue

                    vol += c.volume(ignore_unbounded=ignore_unbounded,
                                    log_coordinate=log_coordinate,
                                    lowerBounds=lowerBounds,
                                    upperBounds=upperBounds,
                                    maxVertices=maxVertices,
                                    limitVertices=limitVertices,
                                    method=method)

                b.extra_columns.children[n].total_volume = vol
            n += 1

    def calculate_total_volumes_geometric(self, cases, b, log_coordinate=True):

        vol_unbounded = 0
        vol_bounded = 0
        for c in cases:
            if c.is_valid() is False or c.is_valid(p_bounds=c.valid_interior_parameter_set()) is False:
                continue
            values = [c.case_number, c.signature]
            values += [self.value_for_extra_column(c, column, filtering_mode=True) for column in b.extra_columns.children]
            if self.show_case(values, b.extra_columns.children, b.filters) is False:
                continue
            vol_unbounded += c.volume_geometric_mean(ignore_unbounded=False, log_coordinate=log_coordinate)
            vol_bounded += c.volume_geometric_mean(ignore_unbounded=True, log_coordinate=log_coordinate)
        vol = [vol_unbounded, vol_bounded]
        return vol

    def save_table(self, b):

        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        controller.save_widget_data(b)

    def generate_multiple_signatures(self, rhs):

        signature = rhs
        controller = self.controller.ds
        i = 0
        siglist = []
        wild_cards = []
        while i < len(signature):
            if signature[i] == '(':
                start = i+1
                while signature[i] != ')':
                    i += 1
                siglist.append(int(signature[start:i]))
            elif signature[i] == '*':
                num_wild = int(controller.signature[i])
                for j in range(num_wild):
                    new_sig = signature.replace('*', str(j+1), 1)
                    wild_cards += self.generate_multiple_signatures(new_sig)
                return wild_cards
            else:
                siglist.append(int(signature[i]))
            i += 1

        return [siglist]

    def show_case_multiple_signatures(self, siglist, a_filter, lhs):
        showCase = True
        for signature_int in siglist:
            rhs = ''.join(str(x) for x in signature_int).replace(" ", "")
            if a_filter[1].value == '==' and lhs != rhs:
                showCase = False
                break
            if a_filter[1].value == '!=' and lhs == rhs:
                showCase = False
                break
        return showCase

    def show_case(self, values, columns, filters):
        showCase = True
        for filter_box in filters.children:
            a_filter = filter_box.children
            column_number = int(a_filter[0].value)
            if column_number == 1:
                rhs = int(a_filter[2].value)
                lhs = int(values[column_number - 1])
            elif column_number == 2:
                rhs = str(a_filter[2].value).replace(" ", "")
                lhs = str(values[column_number - 1]).replace(" ", "")
                if rhs.find("*") != -1:
                    siglist = self.generate_multiple_signatures(rhs)
                    showCase = self.show_case_multiple_signatures(siglist, a_filter, lhs)
                    break
            else:
                column_type = columns[column_number - 3].header.value
                if column_type == 'Validity' or (column_type == '# eigenvalues w/ positive real part' and
                                                 columns[column_number - 3].root_method.value == 'Numpy Functions' and
                                                 columns[column_number - 3].show_complex_conjugates.value is True)\
                        or (column_type == 'Is a Blowup') or (column_type == 'Is False Blowing'):
                    rhs = str(a_filter[2].value)
                    lhs = str(values[column_number - 1])
                else:
                    if values[column_number - 1] == '-' or values[column_number - 1] == '*':
                        showCase = False
                        break
                    rhs = float(a_filter[2].value)
                    lhs = float((str(values[column_number - 1])).replace("*",
                                                                         ""))
            if a_filter[1].value == '==' and lhs != rhs:
                showCase = False
                break
            if a_filter[1].value == '>=' and lhs < rhs:
                showCase = False
                break
            if a_filter[1].value == '<=' and lhs > rhs:
                showCase = False
                break
            if a_filter[1].value == '>' and lhs <= rhs:
                showCase = False
                break
            if a_filter[1].value == '<' and lhs >= rhs:
                showCase = False
                break
            if a_filter[1].value == '!=' and lhs == rhs:
                showCase = False
                break
        return showCase

    def value_for_extra_column(self, case, column, filtering_mode=False):
        controller = self.controller
        s = '<td style="padding:0 15px 0 15px;">'
        if filtering_mode is True:
            volumes = 1
        cyclical = False
        if case.is_cyclical is True:
            cyclical = True
            case = case.original_case
        if case.is_valid() is False:
            if cyclical is True:
                value = '*'
            else:
                value = '-'
            return value
        if column.header.value == 'Validity':
            value = '+'
        elif column.header.value == 'Is a Blowup':
            if case.is_unstable is False:
                value = '-'
            else:
                value = '+'
        elif column.header.value == 'Is False Blowing':
            if case.is_false_blowing is False:
                value = '-'
            else:
                value = '+'
        elif column.header.value == '# eigenvalues w/ positive real part':
            pvals = case.valid_interior_parameter_set()
            if column.root_method.value == 'Routh Array':
                value = case.positive_roots(pvals)
            else:
                    if column.show_complex_conjugates.value is True:
                        if case.has_complex_conjugates(pvals) is True:
                            value = '+'
                        else:
                            value = '-'
                    else:
                        value = case.positive_roots_numpy(pvals)
        elif column.header.value == 'Volume':
            if column.volume_method.value != 'Vertex Enumeration':

                ignore_unbounded = column.ignore_infinity.value
                log_coordinate = column.log_coordinates.value
                method = column.volume_method.value
                lowerBounds = 10 ** column.lower_bound_par_space.value
                upperBounds = 10 ** column.upper_bound_par_space.value

                value = case.volume(ignore_unbounded=ignore_unbounded,
                                    log_coordinate=log_coordinate,
                                    method=method,
                                    lowerBounds=lowerBounds,
                                    upperBounds=upperBounds,
                                    )
            else:
                lb = controller.pvals.copy()
                ub = controller.pvals.copy()
                maxVertices = int(column.nr_vertices.value)
                for key in lb.keys():
                    lb[key] = 10**column.lower_bound_par_space.value    #1E-6
                    ub[key] = 10**column.upper_bound_par_space.value    #1E6
                value, vertices = case.volume_lrs(lb, ub, maxVertices, column.limit_vertices.value)
                value = value if vertices != 0 else 0

        elif column.header.value == 'Dimension':
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value        #1E-6
                ub[key] = 10**column.upper_bound_par_space.value        #1E6
            value = case.dimension(lb, ub)

        elif column.header.value == 'Normalized Volume':
            total_vol = column.total_volume if filtering_mode is False else volumes
            if total_vol == 0:
                total_vol = 1
                print("Warning, setting total volumes to 1")
            value = case.volume(ignore_unbounded=column.ignore_infinity.value,
                                log_coordinate=column.log_coordinates.value,
                                lowerBounds=10**column.lower_bound_par_space.value,
                                upperBounds=10**column.upper_bound_par_space.value,
                                maxVertices=column.nr_vertices.value,
                                limitVertices=column.limit_vertices.value,
                                method=column.volume_method.value,
                                )/total_vol
            value = round_sig(value) if value != 0.0 else 0.0
        else:
            xd = str(column.dependent.value)
            xi = str(column.independent.value)
            value = case.ssystem.log_gain(xd, xi)
        return value

    def html_for_extra_column(self, case, column):
        controller = self.controller
        s = '<td style="padding:0 15px 0 15px;">'
        cyclical = False
        if case.is_cyclical is True:
            cyclical = True
            case = case.original_case
        if case.is_valid() is False:
            if cyclical is True:
                s += '*'
            else:
                s += '-'
            s += '</td>'
            return s
        if column.header.value == 'Validity':
            s += '+'
        elif column.header.value == 'Is a Blowup':
            if case.is_unstable is False:
                s += '-'
            else:
                s += '+'
        elif column.header.value == 'Is False Blowing':
            if case.is_false_blowing is False:
                s += '-'
            else:
                s += '+'
        elif column.header.value == '# eigenvalues w/ positive real part':
            pvals = case.valid_interior_parameter_set()
            if column.root_method.value == 'Routh Array':
                s += str(case.positive_roots(pvals))
            else:
                if column.show_complex_conjugates.value is True:
                        if case.has_complex_conjugates(pvals) is True:
                            s += '+'
                        else:
                            s += '-'
                else:
                        s += str(case.positive_roots_numpy(pvals))

        elif column.header.value == 'Volume':
            if column.volume_method.value != 'Vertex Enumeration':
                ignore_unbounded = column.ignore_infinity.value
                log_coordinate = column.log_coordinates.value
                lowerBounds = 10 ** column.lower_bound_par_space.value
                upperBounds = 10 ** column.upper_bound_par_space.value
                method = column.volume_method.value

                vol = case.volume(ignore_unbounded=ignore_unbounded,
                                  log_coordinate=log_coordinate,
                                  lowerBounds=lowerBounds,
                                  upperBounds=upperBounds,
                                  method=method,
                                  )
                s += '%.2E' % Decimal(str(round_sig(vol))) if vol != 0.0 else '0.0'
            else:
                lb = controller.pvals.copy()
                ub = controller.pvals.copy()
                maxVertices = int(column.nr_vertices.value)
                for key in lb.keys():
                    lb[key] = 10**column.lower_bound_par_space.value       #1E-6
                    ub[key] = 10**column.upper_bound_par_space.value       #1E6
                value, vertices = case.volume_lrs(lb, ub, maxVertices, column.limit_vertices.value)
                value = value if vertices != 0 else 0
                s += '%.2E' % Decimal(str(round_sig(value))) if value != 0 else '0.0'

        elif column.header.value == 'Dimension':
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value           #1E-6
                ub[key] = 10**column.upper_bound_par_space.value           #1E6
            s += str(case.dimension(lb, ub))

        elif column.header.value == 'Normalized Volume':
            total_volume = column.total_volume
            value = case.volume(ignore_unbounded=column.ignore_infinity.value,
                                log_coordinate=column.log_coordinates.value,
                                lowerBounds=10**column.lower_bound_par_space.value,
                                upperBounds=10**column.upper_bound_par_space.value,
                                maxVertices=column.nr_vertices.value,
                                limitVertices=column.limit_vertices.value,
                                method=column.volume_method.value,
                                )/total_volume
            value = round_sig(value) if value != 0.0 else 0.0
            s += '%.2E' % Decimal(str(value))

        else:
            xd = str(column.dependent.value)
            xi = str(column.independent.value)
            s += '%.3f' % case.ssystem.log_gain(xd, xi)
        s += '</td>'
        return s


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)
