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
        case_id = widgets.TextWidget(description='Case', value = '')
        by_signature = widgets.CheckboxWidget(description='By Signature', value=True)
        get_parameters = widgets.ButtonWidget(description='Get Parameters for Case')
        get_parameters.on_click(self.calculate_cases)
        pvals_widgets = [widgets.FloatTextWidget(description=i, value=controller.pvals[i]) 
                         for i in controller.ds.independent_variables]
        get_parameters.case_id = case_id
        get_parameters.by_signature = by_signature
        get_parameters.pvals = pvals_widgets
        wi = widgets.ContainerWidget(children=[case_id, by_signature, get_parameters] + pvals_widgets)
        button = widgets.ButtonWidget(value=False, description='Done')
        button.on_click(self.update_parameters)
        button.pvals = pvals_widgets
        button.wi = wi
        edit_pvals = widgets.ContainerWidget(description='Edit Parameters', children=[wi, button])
        controller.update_child('Edit Parameters', edit_pvals)

    def calculate_cases(self, b):
        controller = self.controller
        case = controller.ds(str(b.case_id.value), by_signature=b.by_signature.value)
        parameters = case.valid_interior_parameter_set()
        if len(parameters) == 0:
            return
        for i in b.pvals:
            i.value = parameters[str(i.description)]
        
    def update_parameters(self, b):
        pvals = b.pvals
        controller = self.controller 
        controller.pvals.update({str(i.description):str(i.value) for i in pvals})
        controller.update_child('Edit Parameters', None)       
