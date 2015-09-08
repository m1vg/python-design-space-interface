import dspace
import dspace.plotutils
import dspace.display

import numpy as np

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


class DisplayTables(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'tables_widget', None)
        
    def create_tables_widget(self):
        
        controller = self.controller
        self.tables_widget = VBox()
        controller.update_child('Tables', self.tables_widget)
        
    def add_table(self, table_data):
        controller = self.controller
        tables = controller.table_data
        new_data = table_data.split('<caption>')
        new_data[1] = 'Table {0}. {1}'.format(len(tables)+1, new_data[1])
        table_data = '<caption>'.join(new_data)
        tables.append(table_data)
        self.add_table_widget(table_data)
        
    def add_table_widget(self, table_data):
        html_widget = HTML(value=table_data)
        wi = Popup(children=[html_widget])
        wi.set_css('height', '300px')
        children = [i for i in self.tables_widget.children]      
        children.append(wi)
        self.tables_widget.children = children
        
    def load_widgets(self):
        controller = self.controller
        for data in controller.table_data:
            self.add_table_widget(data)
                