import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from case_widget import DisplayCase

class CasesTable(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'table', widgets.ContainerWidget())
    
    def cases_table_widget(self):
        controller = self.controller
        options = widgets.RadioButtonsWidget(values = ['None',
                                                       'All',
                                                       'Valid',
                                                       'User Specified'],
                                             value='None')
        cases = widgets.TextareaWidget(description='User specified cases')
        by_signature = widgets.CheckboxWidget(description='Cases indicated by signature?', value=True)
        constraints = widgets.TextareaWidget(description='Biological constraints:')
        add_column = widgets.ButtonWidget(description='Add column')
        remove_column = widgets.ButtonWidget(description='Remove column')
        remove_column.visible = False
        add_column.remove_column = remove_column
        extra_columns = widgets.ContainerWidget(children=[])
        add_column.on_click(self.add_cases_column)
        add_column.column_block = extra_columns
        remove_column.on_click(self.remove_cases_column)
        remove_column.column_block = extra_columns
        case_id = widgets.TextWidget(description='Report for case:')
        wi = widgets.ContainerWidget(children=[options,
                                               cases,
                                               by_signature,
                                               constraints, 
                                               add_column,
                                               extra_columns, 
                                               remove_column
                                               ])
        button = widgets.ButtonWidget(value=False, description='Create/modify cases table')
        button.on_click(self.show_cases)
        button.wi = wi
        button.options = options
        button.cases = cases
        button.by_signature = by_signature
        button.extra_columns = extra_columns
        button.constraints = constraints
        cases_table = widgets.ContainerWidget(description='Cases Table', children=[wi, button, self.table])
        wi.visible = False
        return ('Phenotypic Repertoire', cases_table)
    
    def add_cases_column(self, b):
        controller = self.controller
        children = [i for i in b.column_block.children]
        new_column = [widgets.DropdownWidget(description='Column ' + str(len(children)+3), 
                                               values=['Validity',
                                                       '# eigenvalues w/ positive real part',
                                                       'Log-Gain'],
                                               value='Log-Gain'),
                      widgets.DropdownWidget(description='Log-gain (dependent variable)', 
                                               values=controller.ds.dependent_variables),
                      widgets.DropdownWidget(description='Log-gain (independent variable)', 
                                               values=controller.ds.independent_variables)]
        column = widgets.ContainerWidget(children=new_column)
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
            constraints = [i.strip() for i in constraints]
        else:
            constraints = None
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
            s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td><td style="padding:0 15px 0 15px;">{1}</td>'.format(c.case_number,c.signature)
            for column in b.extra_columns.children:
                s += self.html_for_extra_column(c, column)
            s += '</tr>\n'
        s += '</table><caption>'
        s += 'Note: # of eigenvalues w/ positive real part is calculated using a representative set of parameter values, and may not be reflective of all potential behaviors.'
        s += '</caption></div>'
        html_widget = widgets.HTMLWidget(value = s)
        save_table = widgets.ButtonWidget(description='Save Table')
        save_table.table_data = s
        save_table.on_click(self.save_table)
        table_container = widgets.ContainerWidget(children=[save_table, 
                                                        html_widget])
        table_container.set_css('height', '300px')
        self.table.children = [table_container]
        
    def save_table(self, b):
        
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        
        
    
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