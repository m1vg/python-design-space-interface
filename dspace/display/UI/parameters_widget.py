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
    Popup = HBox
    
from IPython.display import clear_output, display

import matplotlib.pyplot as plt

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  


class EditParameters(object):
        
    def __init__(self, controller):
        setattr(self, 'controller', controller)

    def edit_parameters_widget(self):
        controller = self.controller
        if controller.ds is None:
            return
        pvals_widgets = [FloatText(description=i, value=controller.pvals[i]) 
                         for i in sorted(controller.pvals.keys())]
        wi = VBox(children=pvals_widgets)
        button = Button(value=False, description='Done')
        button.on_click(self.update_parameters)
        button.pvals = pvals_widgets
        button.wi = wi
        edit_pvals = VBox(description='Edit Parameters', children=[wi, button])
        controller.update_child('Edit Parameters', edit_pvals)

        
    def update_parameters(self, b):
        pvals = b.pvals
        controller = self.controller 
        controller.pvals.update({str(i.description):str(i.value) for i in pvals})
        controller.update_child('Edit Parameters', None)       
