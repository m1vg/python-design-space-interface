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
        setattr(self, 'latex', None)
        
    def create_system_widget(self):
        controller = self.controller 
        if controller.ds is None:
            return
        self.html = widgets.HTMLWidget()
        self.latex = widgets.LatexWidget()
        wi = widgets.ContainerWidget(children=[self.html, self.latex])
        self.update_display()
        controller.update_child('System', wi)
        
    def update_display(self):
        controller = self.controller 
        html_string = '<b>Name:</b> ' + ' '.join([controller.name, controller.version]) + '<br>'
        html_string += '<b>Number of Cases:</b> ' + str(controller.ds.number_of_cases) + '<br>'
        html_string += '<b>System Signature:</b> [' + controller.ds.signature + ']<br><hr>'
        html_string += '<b>Equations:</b><br><br>'
        latex_string = controller.ds.equations._repr_latex_()
        self.latex.value = latex_string
        self.html.value = html_string

        
        