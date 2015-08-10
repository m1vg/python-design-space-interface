import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from system_widget import DisplaySystem
from symbols_widget import EditSymbols
from cases_widget import CasesTable
from case_widget import CaseReport
from co_localize_widget import CaseColocalization
from figures_widget import MakePlot, DisplayFigures
from tables_widget import DisplayTables
from parameters_widget import EditParameters

import pickle
import base64

class WidgetSavedData(object):
    
    @staticmethod
    def load_widget_data(interactive):
        
        f = open(interactive.name, 'r')
        saved_data = pickle.load(f)
        f.close()
        figure_data = saved_data.saved['figure_data']
        figure_data = [(base64.b64decode(i[0]), i[1], i[2]) for i in figure_data]
        saved_data.saved['figure_data'] = figure_data
        interactive.__dict__.update(saved_data.saved)
    
    def __init__(self, interactive):
        setattr(self, 'saved', {})
        save_fields = ['ds', 
                       'equations',
                       'name',
                       'pvals',
                       'cyclical',
                       'codominance',
                       'auxiliary',
                       'constraints',
                       'table_data',
                       'symbols',
                       'options',]
        self.saved.update({i:interactive.__dict__[i] for i in save_fields})
        figure_data = interactive.figure_data
        figure_data = [(base64.b64encode(i[0]), i[1], i[2]) for i in figure_data]
        self.saved['figure_data'] = figure_data
        
    def save_data(self):
        
        f = open(self.saved['name'], 'w')
        pickle.dump(self, f)
        f.close() 
        
                       

class InteractiveInput(object):
    
    def __init__(self, name='', equations=None, parameters=None,
                 get_parameters=None, auxiliary_variables=[], constraints=[],
                 symbols={}, resolve_cycles=False, resolve_codominance=False,
                 centered_axes=False, xaxis=None, yaxis=None,
                 x_range=[1e-3, 1e3], y_range=[1e-3, 1e3],
                 zlim=None, by_signature=False, parameter_dict=None,
                 included_cases=None, **kwargs):
        ''' 
        '''
        setattr(self, 'ds', None)
        setattr(self, 'equations', [])
        setattr(self, 'name', name)
        setattr(self, 'pvals', None)
        setattr(self, 'cyclical', resolve_cycles)
        setattr(self, 'codominance', resolve_codominance)
        setattr(self, 'auxiliary', [])
        setattr(self, 'constraints', [])
        setattr(self, 'symbols', symbols)
        setattr(self, 'widget', widgets.TabWidget())
        setattr(self, 'figures', None)
        setattr(self, 'figure_data', [])
        setattr(self, 'tables', None)
        setattr(self, 'table_data', [])
        setattr(self, 'display_system', None)
        setattr(self, 'options', dict(kwargs))
        setattr(self, 'actions', widgets.TabWidget())
        self.options.update(center_axes=centered_axes, 
                            xaxis=xaxis, yaxis=yaxis,
                            range_x=x_range, range_y=y_range, zlim=zlim, 
                            by_signature=by_signature, 
                            parameter_dict=parameter_dict, 
                            get_parameters=get_parameters,
                            included_cases=included_cases)
        if equations is not None:
            self.equations = equations
            if isinstance(equations, list) is False:
                self.equations = [equations]
        
        if auxiliary_variables is not None:
            self.auxiliary = auxiliary_variables
            if isinstance(auxiliary_variables, list) is False:
                self.auxiliary = [auxiliary_variables]
                
        if constraints is not None:
            self.constraints = constraints
            if isinstance(constraints, list) is False:
                self.constraints = [constraints]
                
        if parameters is not None:
            self.pvals = dspace.VariablePool(names = parameters.keys())
            self.pvals.update(parameters)
           
        display(self.widget)   
        self.update_child('Options', self.edit_equations_widget())
        
    def defaults(self, key):
        if key in self.options:
            return self.options[key]
        else:
            return None        
        
    def reload_widgets(self):
        self.widget = widgets.TabWidget()
        display(self.widget)
        self.update_child('Options', self.edit_equations_widget())
        if self.ds is None:
            return
        self.update_widgets()
        self.figures.load_widgets()
        
    def child_with_name(self, name):
        children = self.widget.children
        for i in xrange(len(children)):
            if children[i].description == name:
                return self.widget.children[i]
        return None
    
    def update_child(self, name, child, index=None):
        previous_children = self.widget.children
        children = [i for i in self.widget.children]
        added = False
        for i in xrange(len(children)):
            if children[i].description == name:
                self.widget._titles = {}
                old = children.pop(i)
                break
        if added is False:
            if child is not None:
                child.description = name
                children.append(child)
        self.widget.children = children
        for (i, child) in enumerate(children):
            self.widget.set_title(i, child.description)
    
    def enumerate_phenotypes_menu(self, b):
        cases_table = CasesTable(self)
        return cases_table.cases_table_widget()
        
        
    def case_report_menu(self, b):
        case_report = CaseReport(self, 
                                 self.defaults('by_signature'))
        return case_report.create_case_widget()

    def co_localize_menu(self, b):
        case_report = CaseColocalization(self, 
                                         self.defaults('by_signature'))
        return case_report.colocalization_widget()
    
    def create_plot_menu(self, b):
        plot = MakePlot(self)
        return plot.create_plot_widget()
        
    def load_widget(self, b):
        saved = WidgetSavedData.load_widget_data(self)
        b.wi.visible = True
        if self.ds is None:
            return
        self.make_options_menu(b.button)
        self.update_widgets()
        self.figures.load_widgets()       
        self.tables.load_widgets()       

    def save_widget_data(self, b):
        save = WidgetSavedData(self)
        save.save_data()
        
    def make_options_menu(self, b):
        wi = b.wi
        b.visible = False
        b.name.visible = False
        b.load.visible = False
        actions = [('Enumerate Phenotypes', self.enumerate_phenotypes_menu),
                   ('Analyze Case', self.case_report_menu),
                   ('Co-localize Phenotypes', self.co_localize_menu),
                   ('Create Plot', self.create_plot_menu)]
        actions_h = widgets.HTMLWidget(value='<b>Actions</b>')
        options_h = widgets.HTMLWidget(value='<b>Options</b>')
        options = [('Edit Symbols', self.create_edit_symbols),
                   ('Edit Parameters', self.create_edit_parameters),
                   ('Save widget', self.save_widget_data)]
        actions_w = []
        
        for name, method in actions:
            button = widgets.ButtonWidget(description=name)
            button.on_click(method)
            if method is None:
                button.disabled = True
            actions_w.append(button)
        options_w = []
        for name, method in options:
            button = widgets.ButtonWidget(description=name)
            button.on_click(method)
            if method is None:
                button.disabled = True
            options_w.append(button)
        wi.children = [actions_h, self.actions] + [options_h] + options_w 
        for title, method in actions:
            title, widget = method(self)
            children = [i for i in self.actions.children] + [widget]
            self.actions.children = children
            self.actions.set_title(len(children)-1, title)
        ##     
        ## self.enumerate_phenotypes_menu(self)
        ## self.case_report_menu(self)
        ## self.co_localize_menu(self)
        ## self.create_plot_menu(self)
        
        
    def update_widgets(self):
        self.display_system = DisplaySystem(self)
        self.display_system.create_system_widget()
        self.figures = DisplayFigures(self)
        self.figures.create_figures_widget()
        self.tables = DisplayTables(self)
        self.tables.create_tables_widget()
                
    def edit_equations_widget(self):
        parameter_dict = self.options['parameter_dict']
        if parameter_dict is None:
            parameter_dict = {} 
        name = widgets.TextWidget(description='* Name', value=self.name)
        equations=widgets.TextareaWidget(description='* Equations',
                                         value='\n'.join(self.equations))
        aux=widgets.TextareaWidget(description='Auxiliary Variables',
                                   value=', '.join(self.auxiliary))
        html = widgets.HTMLWidget(value='<b>Architectural Constraints</b>')
        constraints=widgets.TextareaWidget(description='Parameters',
                                           value=', '.join(self.constraints)
                                           )
        options_html = widgets.HTMLWidget(value='<b>Additional Options</b>')
        cyclical = widgets.CheckboxWidget(description='Check for Cycles',
                                          value = self.cyclical)
        codominance = widgets.CheckboxWidget(description='Check for Co-dominance',
                                             value = self.codominance)
        replacements=widgets.TextareaWidget(description='Kinetic Orders',
                                            value=', '.join(
                                             [key + '=' + value for (key,value) in
                                             parameter_dict.iteritems()]))
        wi = widgets.ContainerWidget(children=[equations, 
                                               aux, html,
                                               constraints, replacements,
                                               options_html, cyclical,
                                               codominance,
                                               ])
        button = widgets.ButtonWidget(value=False, description='Create Design Space')
        button.on_click(self.make_design_space)
        button.equations = equations
        button.aux = aux
        button.constraints = constraints
        button.cyclical = cyclical
        button.codominance = codominance
        button.replacements = replacements
        button.wi = wi
        button.name = name
        load = widgets.ButtonWidget(value=False, description='Load Widget')
        load.on_click(self.load_widget)
        load.equations = equations
        load.aux = aux
        load.constraints = constraints
        load.cyclical = cyclical
        load.codominance = codominance
        load.replacements = replacements
        load.wi = wi
        load.name = name
        load.button = button
        button.load = load
        edit_symbols = widgets.ButtonWidget(value=False, description='Edit Symbols')
        edit_symbols.on_click(self.create_edit_symbols)
        edit_parameters = widgets.ButtonWidget(value=False, description='Edit Parameters')
        edit_parameters.on_click(self.create_edit_parameters)
        button.edit_symbols = edit_symbols
        button.edit_parameters = edit_parameters
        edit_equations = widgets.ContainerWidget(description='Edit Equations', 
                                                 children=[name, wi, 
                                                           edit_symbols,
                                                           edit_parameters,
                                                           button,
                                                           load])
        wi.visible = False
        edit_symbols.visible = False
        edit_parameters.visible = False
        if self.ds is not None:
            wi.visible = True
            self.update_widgets()
            self.make_options_menu(button)
        return edit_equations  
        
    def make_design_space(self, b):
        if b.wi.visible == False:
            b.wi.visible = True
            b.description = 'Done'
            return
        self.equations = [i.strip() for i in str(b.equations.value).split('\n') if len(i.strip()) > 0]
        self.auxiliary = [i.strip() for i in str(b.aux.value).split(',') if len(i.strip()) > 0] 
        self.constraints = [i.strip() for i in str(b.constraints.value).split(',') if len(i.strip()) > 0] 
        self.cyclical = b.cyclical.value
        self.codominance = b.codominance.value
        eq = dspace.Equations(self.equations,
                              auxiliary_variables=self.auxiliary, 
                              latex_symbols=self.symbols)
        self.name = b.name.value
        constraints = self.constraints
        replacements = str(b.replacements.value).strip()
        if len(replacements) > 0:
            replacements = str(b.replacements.value).split(',')
            replacements = [[j.strip() for j in i.split('=')] for i in replacements] 
            parameter_dict = {i:j for i,j in replacements}
        else:
            parameter_dict = {}
        if len(constraints) == 0:
            constraints = None
        self.ds = dspace.DesignSpace(eq, name=b.name.value, 
                                     constraints=constraints,
                                     resolve_cycles=self.cyclical, 
                                     resolve_codominance=self.codominance,
                                     parameter_dict=parameter_dict)
        if self.pvals == None:
            self.pvals = dspace.VariablePool(
                          names=self.ds.independent_variables
                          )
        else:
            pvals = dspace.VariablePool(
                     names=self.ds.independent_variables
                     )
            for i in pvals:
                if i in self.pvals:
                    pvals[i] = self.pvals[i]
            self.pvals = pvals
        get_parameters = self.defaults('get_parameters')
        if get_parameters is not None:
             case = self.ds(get_parameters)
             self.pvals = case.valid_interior_parameter_set()
        self.symbols.update({i:i for i in self.ds.dependent_variables + self.ds.independent_variables if i not in self.symbols})
        self.ds.update_latex_symbols(self.symbols)
        self.update_widgets()
        self.make_options_menu(b)
    
    def create_edit_symbols(self, b):
        Edit = EditSymbols(self)
        Edit.edit_symbols_widget() 
               
    def create_edit_parameters(self, b):
        Edit = EditParameters(self)
        Edit.edit_parameters_widget()
        
