import dspace
import dspace.plotutils
import dspace.display
import pandas as pd

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
        pvals_widgets = [FloatText(description=i,
                                   value=controller.pvals[i])
                         for i in sorted(controller.pvals.keys())]

        identity_widgets = [Dropdown(description=i,
                                    value='K',

                                    values=['K',
                                            'b',
                                            'a',
                                            'a_max',
                                            'a_min',
                                            'Other'],

                                    options=['K',
                                             'b',
                                             'a',
                                             'a_max',
                                             'a_min',
                                             'Other'])
                            for i in sorted(controller.pvals.keys())]

        if controller.pidentity is not None:
            self.update_identity_widgets(identity_widgets)

        tx = HTML(value='For population dynamics studies, please proof the identity of each parameter below. '
                        'Five different identities are possible: K, b (&#946), a_max (&#945_max), '
                        'a_min (&#945_min) and "other".')

        ch = pvals_widgets + [tx] + identity_widgets

        wi = VBox(children=ch)
        button = Button(value=False, description='Done')
        button.on_click(self.update_parameters)
        button.pvals = pvals_widgets
        button.pidentity = identity_widgets
        button.wi = wi

        load_data = Button(value=False, description='Get Values From .xlsx File')
        load_data.controller = controller
        load_data.pvals = pvals_widgets
        load_data.on_click(self.load_parameters_from_file)
        edit_pvals = VBox(description='Edit Parameters', children=[wi, load_data, button])
        controller.update_child('Edit Parameters', edit_pvals)
        
    def update_parameters(self, b):
        pvals = b.pvals
        controller = self.controller 
        controller.pvals.update({str(i.description):str(i.value) for i in pvals})
        self.update_parameter_identity(b)
        controller.update_child('Edit Parameters', None)

    def load_parameters_from_file(self, b):
        # This function should:
        # 1. Generate dataFrame from excel file
        # 2. Loop over each row of dataFrame and update value, if present of field .pvals_widgets

        controller = b.controller
        df = pd.read_excel(controller.name + '_parameters.xlsx')
        par_keys_excel = df[df.columns[0]].values
        par_keys_ds = sorted(controller.pvals.keys())
        header0 = df.columns[0]
        header1 = df.columns[1]
        for key_excel in par_keys_excel:
            if key_excel in par_keys_ds:
                index = par_keys_ds.index(key_excel)
                df_row = df[df[header0].values == key_excel]
                b.pvals[index].value = df_row[header1].values[0]

    def update_parameter_identity(self, b):
        identity = b.pidentity
        controller = self.controller
        controller.pidentity = controller.pvals.copy()
        for i in identity:
            if i.value == 'K':
                value = 1
            elif i.value == 'b':
                value = 2
            elif i.value == 'a_max':
                value = 3
            elif i.value == 'a_min':
                value = 4
            elif i.value == 'a':
                value = 5
            else:
                value = 6
            controller.pidentity.update({str(i.description): value})

    def update_identity_widgets(self, identity_widgets):
        identity = self.controller.pidentity
        for i in identity_widgets:

            parameter = i.description
            numerical_key = identity[parameter]

            if numerical_key == 1:
                par_s = 'K'
            elif numerical_key == 2:
                par_s = 'b'
            elif numerical_key == 3:
                par_s = 'a_max'
            elif numerical_key == 4:
                par_s = 'a_min'
            elif numerical_key == 5:
                par_s = 'a'
            else:
                par_s = 'Other'
            i.value = par_s


