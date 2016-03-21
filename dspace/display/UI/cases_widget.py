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
        options = Select(values = ['None',
                                   'All',
                                   'Valid',
                                   'User Specified'],
                         options = ['None',
                                    'All',
                                    'Valid',
                                    'User Specified'],
                         value='None')
        if 'biological_constraints' not in controller.options:
            bio_constraints = ''
        else:
            bio_constraints = ', '.join(controller.defaults('biological_constraints')) 
        cases = Textarea(description='User specified cases')
        by_signature = Checkbox(description='Cases indicated by signature?', value=True)
        constraints = Textarea(description='Biological constraints:',
                                             value=bio_constraints)
        add_column = Button(description='Add column')
        add_filter = Button(description='Add filter')
        remove_column = Button(description='Remove column')
        remove_column.visible = False
        add_column.remove_column = remove_column
        extra_columns = VBox(children=[])
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
                            add_column,
                            extra_columns, 
                            remove_column,
                            filters,
                            add_filter
                           ])
        button = Button(value=False, description='Create/modify cases table')
        button.on_click(self.show_cases)
        button.wi = wi
        button.options = options
        button.cases = cases
        button.by_signature = by_signature
        button.extra_columns = extra_columns
        button.constraints = constraints
        button.filters = filters
        cases_table = VBox(description='Cases Table', children=[wi, button, self.table])
        wi.visible = False
        return ('Phenotypic Repertoire', cases_table)
    
    def add_cases_filter(self, b):
        controller = self.controller
        children = [i for i in b.filters.children]
        new_filter = [Dropdown(description='Column #', 
                               values=[str(i+1) for i in range(2+len(b.column_block.children))],
                               options=[str(i+1) for i in range(2+len(b.column_block.children))],
                               value='1'),
                      Dropdown(description='Condition', 
                               values=['==', '>', '<', '>=', '<='],
                               options=['==', '>', '<', '>=', '<=']),
                      Text(description='Value'),
                      Button(description='X')]
        new_filter[3].on_click(self.remove_cases_filter)
        current = HBox(children=new_filter)
        new_filter[3].current = current
        new_filter[3].filters = b.filters
        children.append(current)
        b.filters.children = children
    
    #def update_cases_filter(self, filters, columns):
        
    def remove_cases_filter(self, b):
        controller = self.controller
        if len(b.filters.children) == 0:
            return
        children = [i for i in b.filters.children]
        for i,child in enumerate(children):
            if b.current == child:
                children.pop(i)
        b.filters.children = children
        
    def add_cases_column(self, b):
        controller = self.controller
        children = [i for i in b.column_block.children]
        new_column = [Dropdown(description='Column ' + str(len(children)+3), 
                               values=['Validity',
                                       '# eigenvalues w/ positive real part',
                                       'Log-Gain'],
                               options=['Validity',
                                        '# eigenvalues w/ positive real part',
                                        'Log-Gain'],
                               value='Log-Gain'),
                      Dropdown(description='Log-gain (dependent variable)', 
                               values=controller.ds.dependent_variables,
                               options=controller.ds.dependent_variables,),
                      Dropdown(description='Log-gain (independent variable)', 
                               values=controller.ds.independent_variables,
                               options=controller.ds.independent_variables)]
        column = VBox(children=new_column)
        column.header = new_column[0]
        column.dependent = new_column[1]
        column.independent = new_column[2]
        children.append(column)
        b.column_block.children = children
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
            b.description = 'Create/modify cases table'
            return
        if b.wi.visible == False:
            self.table.children = []
            b.wi.visible = True
            b.description = 'Done'
            return 
        self.create_case_table(b, mode=str(b.options.value))
        b.wi.visible = False
        b.description = 'Create/modify cases table'
    
    def create_case_table(self, b, mode='All'):
        controller = self.controller
        constraints = str(b.constraints.value)
        if constraints != '':
            constraints = constraints.split(',')
            constraints = [i.strip() for i in constraints if len(i.strip()) > 0]
        else:
            constraints = []
        controller.set_defaults('biological_constraints', constraints)
        if mode == 'None':
            self.table.children=[]
            return
        s = '<div><table>\n<caption>Cases in the system design space. </caption>\n'
        s += '<tr align=center><td style="padding:0 15px 0 15px;"><b>{0}</b></td><td style="padding:0 15px 0 15px;"><b>{1}</b></td>'.format('  Case Number  ', '  Case Signature  ')
        for column in b.extra_columns.children:
            header = str(column.header.value)
            if header == 'Log-Gain':
                header = 'L('+str(column.dependent.value)+','+str(column.independent.value)+')'
            s += '<td><b>' + header + '</b></td>'
        s += '</tr>'
        if mode == 'Valid':    
            cases = controller.ds(controller.ds.valid_cases(),
                                  constraints=constraints)
        elif mode == 'All':
            cases = controller.ds(
                     [i for i in xrange(1, controller.ds.number_of_cases+1)],
                     constraints=constraints)
        else:
            case_signatures = [i.strip() for i in str(b.cases.value).split(',') if len(i.strip()) > 0] 
            cases = controller.ds([i for i in case_signatures], 
                                  by_signature=b.by_signature.value,
                                  constraints=constraints)
            
        for c in cases:
            #self.show_case(c, b.extra_columns)
            values = [c.case_number, c.signature]
            values += [self.value_for_extra_column(c, column) for column in b.extra_columns.children]
            if self.show_case(values, b.extra_columns.children, b.filters) is False:
                continue
            s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td><td style="padding:0 15px 0 15px;">{1}</td>'.format(c.case_number,c.signature)
            for column in b.extra_columns.children:
                s += self.html_for_extra_column(c, column)    
            s += '</tr>\n'
            
        s += '</table><caption>'
        s += 'Note: # of eigenvalues w/ positive real part is calculated using a representative set of parameter values, and may not be reflective of all potential behaviors.'
        s += '</caption></div>'
        html_widget = HTML(value = s)
        save_table = Button(description='Save Table')
        save_table.table_data = s
        save_table.on_click(self.save_table)
        table_container = VBox(children=[save_table, 
                                         html_widget])
        if ipy_old is True:
            table_container.set_css('height', '300px')
        else:
            table_container.height = '300px'
            table_container.overflow_x = 'auto'
            table_container.overflow_y = 'auto'
        self.table.children = [table_container]
        
    def save_table(self, b):
        
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        controller.save_widget_data(b)
        
    def show_case(self, values, columns, filters):
        showCase = True

        for filter_box in filters.children:
            a_filter = filter_box.children
            column_number = int(a_filter[0].value)
            if column_number == 1:
                rhs = int(a_filter[2].value)
                lhs = int(values[column_number-1])
            elif column_number == 2:
                rhs = str(a_filter[2].value)
                lhs = str(values[column_number-1])
            else:
                column_type = columns[column_number-3].header.value
                if column_type == 'Validity':
                    rhs = str(a_filter[2].value)
                    lhs = str(values[column_number-1])
                else:
                    if values[column_number-1] == '-' or values[column_number-1] == '*':
                        showCase = False
                        break
                    rhs = float(a_filter[2].value)
                    lhs = float(values[column_number-1])
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
        return showCase
    
    def value_for_extra_column(self, case, column):
        controller = self.controller
        s = '<td style="padding:0 15px 0 15px;">'
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
        elif column.header.value == '# eigenvalues w/ positive real part':
            value = case.positive_roots(controller.pvals)
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
        elif column.header.value == '# eigenvalues w/ positive real part':
            s += str(case.positive_roots(controller.pvals))
        else:
            xd = str(column.dependent.value)
            xi = str(column.independent.value)
            s += '%.3f' % case.ssystem.log_gain(xd, xi)
        s += '</td>'
        return s