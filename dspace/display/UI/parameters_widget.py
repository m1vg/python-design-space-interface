import dspace
import dspace.plotutils
import dspace.display

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

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
        pvals_widgets = [widgets.FloatTextWidget(description=i, value=controller.pvals[i]) 
                         for i in sorted(controller.pvals.keys())]
        wi = widgets.ContainerWidget(children=pvals_widgets)
        button = widgets.ButtonWidget(value=False, description='Done')
        button.on_click(self.update_parameters)
        button.pvals = pvals_widgets
        button.wi = wi
        edit_pvals = widgets.ContainerWidget(description='Edit Parameters', children=[wi, button])
        controller.update_child('Edit Parameters', edit_pvals)

        
    def update_parameters(self, b):
        pvals = b.pvals
        controller = self.controller 
        controller.pvals.update({str(i.description):str(i.value) for i in pvals})
        controller.update_child('Edit Parameters', None)       
