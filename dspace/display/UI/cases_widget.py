import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  


class CasesTable(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'table', widgets.ContainerWidget())
    
    def cases_table_widget(self):
        controller = self.controller
#         html_widget = widgets.HTMLWidget(value='')
#         table = widgets.ContainerWidget(children=[html_widget])
#         table.visible = False
        options = widgets.RadioButtonsWidget(values = ['None',
                                                       'All',
                                                       'Valid',
                                                       'User Specified'],
                                             value='None')
        cases = widgets.TextareaWidget(description='User specified cases')
        by_signature = widgets.CheckboxWidget(description='Specify cases by signature', value=True)
        add_column = widgets.ButtonWidget(description='Add column')
        extra_columns = widgets.ContainerWidget(children=[])
        add_column.on_click(self.add_cases_column)
        add_column.column_block = extra_columns
        cases.placeholder='Insert cases by signature, seperated by commas, for user specified cases.  Wildcards are allowed.'
        wi = widgets.ContainerWidget(children=[options, cases, by_signature, add_column, extra_columns])
        button = widgets.ButtonWidget(value=False, description='Show Cases Options')
        button.on_click(self.show_cases)
#         button.html = html_widget
        button.wi = wi
#         button.table = table
        button.options = options
        button.cases = cases
        button.by_signature = by_signature
        button.extra_columns = extra_columns
        cases_table = widgets.ContainerWidget(description='Cases Table', children=[wi, button, self.table])
        controller.update_child('Cases', cases_table)
        wi.visible = False
    
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
        
    def show_cases(self, b):
        controller = self.controller
        if controller.ds == None:
            b.wi.visible = False
            self.table.children = []
            b.description = 'Show Cases Options'
            return
        if b.wi.visible == False:
            self.table.children = []
            b.wi.visible = True
            b.description = 'Done'
            return 
        self.create_case_table(b, mode=str(b.options.value))
        b.wi.visible = False
        b.description = 'Show Cases Options'
    
    def create_case_table(self, b, mode='All'):
        controller = self.controller
        if mode == 'None':
            self.table.children=[]
            return
        s = '<div><table>\n<caption>Table 1. Cases in the system design space. </caption>\n'
        s += '<tr align=center><td style="padding:0 15px 0 15px;"><b>{0}</b></td><td style="padding:0 15px 0 15px;"><b>{1}</b></td>'.format('  Case Number  ', '  Case Signature  ')
        for column in b.extra_columns.children:
            header = str(column.header.value)
            if header == 'Log-Gain':
                header = 'L('+str(column.dependent.value)+','+str(column.independent.value)+')'
            s += '<td><b>' + header + '</b></td>'
        s += '</tr>'
        if mode == 'Valid':    
            cases = controller.ds(controller.ds.valid_cases())
        elif mode == 'All':
            cases = controller.ds([i for i in xrange(1, controller.ds.number_of_cases+1)])
        else:
            case_signatures = [i.strip() for i in str(b.cases.value).split(',') if len(i.strip()) > 0] 
            cases = controller.ds([i for i in case_signatures], by_signature=b.by_signature.value)
        for c in cases:
            s += '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td><td style="padding:0 15px 0 15px;">{1}</td>'.format(c.case_number,c.signature)
            for column in b.extra_columns.children:
                s += self.html_for_extra_column(c, column)
            s += '</tr>\n'
        s += '</table></div>'
        html_widget = widgets.HTMLWidget(value = s)
        table_container = widgets.PopupWidget(children=[html_widget])
        table_container.set_css('height', '300px')
        self.table.children = [table_container]
        
        
    
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