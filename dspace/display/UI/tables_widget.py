import dspace
import dspace.plotutils
import dspace.display

import numpy as np
from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML, Latex

class DisplayTables(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'tables_widget', None)
        
    def create_tables_widget(self):
        
        controller = self.controller
        self.tables_widget = widgets.ContainerWidget()
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
        html_widget = widgets.HTMLWidget(value=table_data)
        wi = widgets.PopupWidget(children=[html_widget])
        wi.set_css('height', '300px')
        children = [i for i in self.tables_widget.children]      
        children.append(wi)
        self.tables_widget.children = children
        
    def load_widgets(self):
        controller = self.controller
        for data in controller.table_data:
            self.add_table_widget(data)
                