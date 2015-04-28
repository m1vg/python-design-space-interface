import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  

from symbols_widget import EditSymbols
from cases_widget import CasesTable
from figures_widget import MakePlot, DisplayFigures
from parameters_widget import EditParameters

class InteractiveInput(object):
    
    def __init__(self, name='', equations=None, parameters=None, auxiliary_variables=[], constraints=[], symbols={}, 
                 resolve_cycles=False, resolve_codominance=False, **kwargs):
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
        children = [i for i in self.widget.children]
        added = False
        for i in xrange(len(children)):
            if children[i].description == name:
                added = True
                children.pop(i)
                if child is not None:
                    child.description = name
                    children.insert(i, child)
                    self.widget.set_title(i, name)
        if added is False:
            if child is not None:
                child.description = name
                children.append(child)
                self.widget.set_title(len(children)-1, name)
        self.widget.children = children
        
    def update_widgets(self):
        cases_table = CasesTable(self)
        cases_table.cases_table_widget()
        plot = MakePlot(self)
        plot.create_plot_widget()
        self.figures = DisplayFigures(self)
        self.figures.create_figures_widget()

        
    def edit_equations_widget(self):
        name = widgets.TextWidget(description='Name', value=self.name)
        equations=widgets.TextareaWidget(description='Equations',
                                         value='\n'.join(self.equations))
        equations.placeholder='Insert equations, one per line.'
        aux=widgets.TextareaWidget(description='Auxiliary Variables',
                                                   value=', '.join(self.auxiliary))
        aux.placeholder='Insert names of auxiliary variables, if any, seperated by commas.'
        constraints=widgets.TextareaWidget(description='Constraints',
                                           value=', '.join(self.constraints))
        constraints.placeholder='Insert constraint inequalities, if any, seperated by commas.'
        cyclical = widgets.CheckboxWidget(description='Check for Cycles',
                                          value = self.cyclical)
        codominance = widgets.CheckboxWidget(description='Check for Co-dominance',
                                             value = self.codominance)
        wi = widgets.ContainerWidget(children=[name, equations,
                                               aux, constraints, cyclical, codominance])
        if self.ds is None:
            description = 'Create Design Space'
        else:
            description = 'Edit Design Space'
        button = widgets.ButtonWidget(value=False, description=description)
        button.on_click(self.make_design_space)
        button.equations = equations
        button.aux = aux
        button.constraints = constraints
        button.cyclical = cyclical
        button.codominance = codominance
        button.wi = wi
        button.name = name
        edit_symbols = widgets.ButtonWidget(value=False, description='Edit Symbols')
        edit_symbols.on_click(self.create_edit_symbols)
        edit_parameters = widgets.ButtonWidget(value=False, description='Edit Parameters')
        edit_parameters.on_click(self.create_edit_parameters)
        button.edit_symbols = edit_symbols
        button.edit_parameters = edit_parameters
        edit_equations = widgets.ContainerWidget(description='Edit Equations', 
                                                 children=[wi, 
                                                           edit_symbols,
                                                           edit_parameters,
                                                           button])
        wi.visible = False
        edit_symbols.visible = False
        edit_parameters.visible = False
        if self.ds is not None:
            edit_symbols.visible = True
            edit_parameters.visible = True    
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
        if len(constraints) == 0:
            constraints = None
        self.ds = dspace.DesignSpace(eq, name=b.name.value, constraints=constraints,
                                     resolve_cycles=self.cyclical, 
                                     resolve_codominance=self.codominance)
        if self.pvals == None:
            self.pvals = dspace.VariablePool(names=self.ds.independent_variables)
        else:
            pvals = dspace.VariablePool(names=self.ds.independent_variables)
            for i in pvals:
                if i in self.pvals:
                    pvals[i] = self.pvals[i]
            self.pvals = pvals
        self.update_widgets()
        b.wi.visible = False
        b.description = 'Edit Design Space'
        b.edit_symbols.visible = True
        b.edit_parameters.visible = True
    
    def create_edit_symbols(self, b):
        Edit = EditSymbols(self)
        Edit.edit_symbols_widget()
        
    def create_edit_parameters(self, b):
        Edit = EditParameters(self)
        Edit.edit_parameters_widget()
        
