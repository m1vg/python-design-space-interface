import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

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
        self.html = widgets.HTMLWidget()
        self.html_equations = widgets.HTMLWidget()
        self.latex = widgets.LatexWidget()
        self.checkbox = widgets.CheckboxWidget(description='Typeset Equations?',
                                               value=True)
        self.checkbox.on_trait_change(self.changed_check_box, 'value')
        wi = widgets.ContainerWidget(children=[self.html,
                                               self.checkbox,
                                               self.html_equations,
                                               self.latex])
        self.update_display()
        controller.update_child('System', wi)
        
    def update_display(self):
        controller = self.controller 
        html_string = '<b>Name:</b> ' + ' '.join([controller.name, controller.version]) + '<br>'
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
        