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


class DisplaySystem(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'html', None)
        setattr(self, 'html_equations', None)
        setattr(self, 'latex', None)

    def create_system_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        self.html = HTML()
        self.html_equations = HTML()
        self.latex = Latex()
        self.checkbox = Checkbox(description='Typeset Equations?',
                                 value=True)
        self.checkbox.on_trait_change(self.changed_check_box, 'value')
        wi = VBox(children=[self.html,
                            self.checkbox,
                            self.html_equations,
                            self.latex,
                            ])
        self.update_display()
        controller.update_child('System', wi)
        
    def update_display(self):
        controller = self.controller 
        html_string = '<b>Name:</b> ' + controller.name
        html_string += ' (<font color="darkred">v{0}</font>)<br>'.format(controller.version) if controller.version != '' else '<br>'
        html_string += '<b>Number of Cases:</b> ' + str(controller.ds.number_of_cases) + '<br>'
        html_string += '<b>System Signature:</b> [' + controller.ds.signature + ']<br><hr>'
        html_string += '<b>Equations:</b><br><br>'
        if self.checkbox.value is True:
            self.html_equations.value = ''
            latex_string = controller.ds.equations._repr_latex_()
        else:
            self.html_equations.value = '<br>'.join([str(i) for i in controller.ds.equations]) 
            latex_string = ''
        self.latex.value = latex_string
        self.html.value = html_string
    
    def changed_check_box(self, name, value):
        self.update_display()
        