import dspace
import dspace.plotutils
import dspace.display
import numpy as np
import cPickle as pickle

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
    ipy_old = True
else:
    from ipywidgets import *
    from popup import Popup
    ipy_old = False

    
from IPython.display import clear_output, display


class DisplayTables(object):
    
    def __init__(self, controller):
        setattr(self, 'controller', controller)
        setattr(self, 'tables_widget', None)

        
    def create_tables_widget(self):
        
        controller = self.controller
        self.tables_widget = VBox()
        self.delete_tables_button = Button(description="Delete Selected Tables")
        self.delete_tables_button.on_click(self.open_delete_menu)
        self.delete_menu = VBox(children=[self.delete_tables_button])

        self.tables = VBox()
        self.tables.children = [self.delete_menu, self.tables_widget]

        controller.update_child('Tables', self.tables)
        
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
        if ipy_old is True:
            wi.set_css('height', '300px')
        children = [i for i in self.tables_widget.children]      
        children.append(wi)
        self.tables_widget.children = children
        
    def load_widgets(self):
        controller = self.controller
        for data in controller.table_data:
            self.add_table_widget(data)

    def open_delete_menu(self, b):  # Miguel
        description = Latex()
        description.value = "Specify tables to delete. Please note that deleted items cannot be recovered!"
        indText = Textarea(description="Individual tables (comma separated values)")
        rangeText = Textarea(description="Range of tables (e.g 1-2)")
        deleteButton = Button(description="Delete tables")
        cancelButton = Button(description="Cancel")
        self.delete_menu.children = [description, indText, rangeText, deleteButton, cancelButton]
        deleteButton.indValue = indText
        deleteButton.rangeValue = rangeText
        deleteButton.on_click(self.updatefile)
        cancelButton.on_click(self.cancelDelete)

    def cancelDelete(self, b):
        self.delete_menu.children = [self.delete_tables_button]

    def updatefile(self, b):  # Miguel

        # load file dsipy
        controller = self.controller
        version = controller.version
        if version == '':
            file_name = controller.name + ".dsipy"
        else:
            file_name = controller.name + '-V' + str(version) + ".dsipy"
        f = open(file_name, "r")
        saved_data = pickle.load(f)
        f.close()

        tables = saved_data.saved["table_data"]

        # Extract Single Figures
        string = b.indValue.value
        string = string.split(',')
        if b.indValue.value:
            tabID = [int(s) for s in string]
        else:
            tabID = []

        if len(tabID) > 0:
            for s in tabID:
                del tables[s - 1]
                del controller.table_data[s - 1]

        # Extract Range
        string = b.rangeValue.value
        string = string.split('-')
        if b.rangeValue.value:
            indicesRange = [int(s) for s in string]
        else:
            indicesRange = []

        if len(indicesRange) == 2:
            del tables[indicesRange[0] - 1:indicesRange[1]]
            del controller.table_data[indicesRange[0] - 1:indicesRange[1]]

        ####  update tables list:
        tables_updated_nr = []
        if len(tables) > 0:
            count = 1
            for row_i in tables:
                old_table_nr = row_i.split('<div><table>\n<caption>')[1].split('.')[0]
                new_table_nr = 'Table ' + str(count)
                updated_row_i = row_i.replace(old_table_nr, new_table_nr, 1)
                tables_updated_nr.append(updated_row_i)
                count += 1
        tables = tables_updated_nr

        saved_data.saved["table_data"] = tables
        controller.table_data = tables
        # controller.figure_data = Figures

        f = open(file_name, "w")
        pickle.dump(saved_data, f)
        f.close()
        sucessMessage = Latex()
        sucessMessage.value = "Tables were deleted and file was updated!"
        self.delete_menu.children = [sucessMessage, self.delete_tables_button]

        controller.tables.create_tables_widget()
        self.load_widgets()
                