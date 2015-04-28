import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg 


class EditSymbols(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        
    def edit_symbols_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        symbols = {i:controller.symbols[i] for i in controller.symbols if i in controller.ds.dependent_variables + controller.ds.independent_variables}
        symbols.update({i:i for i in controller.ds.dependent_variables + controller.ds.independent_variables if i not in symbols})
        self.symbols = symbols
        symbols_widgets = [widgets.TextWidget(description=i, value=self.symbols[i]) for i in symbols]
        wi = widgets.ContainerWidget(children=symbols_widgets)
        button = widgets.ButtonWidget(value=False, description='Done')
        button.on_click(self.update_symbols)
        button.symbols = symbols_widgets
        button.wi = wi
        edit_symbols = widgets.ContainerWidget(description='Edit Symbols', children=[wi, button])
        controller.update_child('Edit Symbols', edit_symbols)
        wi.visible = True
    
    def update_symbols(self, b):
        symbols = b.symbols
        controller = self.controller 
        if b.wi.visible == False:
            for i in symbols:
                i.value = self.symbols[str(i.description)]
            b.wi.visible = True
            b.description = 'Done'
            return 
        controller.symbols.update({str(i.description):str(i.value) for i in symbols})
        controller.ds.update_latex_symbols(self.symbols)
        controller.update_child('Edit Symbols', None)
        controller.display_system.update_display()