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
    from IPython.html.widgets import FloatTextWidget as FloatText
    from IPython.html.widgets import ImageWidget as Image
    VBox = Box
    HBox = Box
else:
    from ipywidgets import *
    
from IPython.display import clear_output, display

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
        symbols_widgets = [Text(description=i, value=self.symbols[i]) for i in symbols]
        wi = VBox(children=symbols_widgets)
        button = Button(value=False, description='Done')
        button.on_click(self.update_symbols)
        button.symbols = symbols_widgets
        button.wi = wi
        edit_symbols = VBox(description='Edit Symbols', children=[wi, button])
        controller.update_child('Edit Symbols', edit_symbols)
        wi.visible = True
    
    def update_symbols(self, b):
        symbols = b.symbols
        controller = self.controller 
        controller.symbols.update({str(i.description):str(i.value) for i in symbols})
        controller.ds.update_latex_symbols(controller.symbols)
        controller.update_child('Edit Symbols', None)
        controller.display_system.update_display()
        
        