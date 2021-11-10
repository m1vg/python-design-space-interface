from __future__ import division

import itertools

import numpy as np
import matplotlib as mt
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from math import log10, floor, ceil

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.expressions import Expression
from dspace.models.cyclicalcase import CyclicalCase
import dspace.models.case
import dspace.models.designspace

import dspace.plotutils.designspace_plot

from dspace.plotutils.monkey_patching import monkeypatch_method
from dspace.models.designspace import DesignSpace

import dspace.plotutils.case_plot
from dspace.models.designspace import sort_cases

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
    from IPython.html.widgets import FloatSliderWidget as FloatSlider
    from IPython.html.widgets import ImageWidget as Image
    VBox = Box
    HBox = Box
else:
    from ipywidgets import *
    Popup = HBox
    
from IPython.display import clear_output, display
import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg


def colorbar_tabs_html(colors, height=None):
    tab_dicts = {}
    html_widgets = []
    for i in colors:
        key = len(i.split(','))
        if key not in tab_dicts:
            tab_dicts[key] = {}
        tab_dicts[key][i] = '#%02x%02x%02x' % tuple([j * 255 for j in colors[i][:3]])
    keys = sorted(tab_dicts)
    labels = [sorted(tab_dicts[i]) for i in keys]
    lengths = [len(tab_dicts[i]) for i in keys]
    max_length = max(lengths)
    html_str = '<table style="border:0;"' + ('height="' + height + '">' if height is not None else '>')
    for i in xrange(max_length):
        html_str += '<tr style="border:0;">'
        for j in xrange(len(labels)):
            if i < lengths[j]:
                key = keys[j]
                label = labels[j][i]
                html_str += '<td style="border:0; background-color:{0}; padding:10px;" />'.format(tab_dicts[key][label])
                html_str += '<td style="border:0; white-space: nowrap; font-size:100%">' + label + '</td>'
            else:
                html_str += '<td style="border:0;" />'
                html_str += '<td style="border:0;" />'
        html_str += '</tr>'
    html_str += '</table>'
    return html_str


def colorbar_tabs(colors):
    html_str = colorbar_tabs_html(colors)
    return html_str

def make_2D_slice(ds=None, p_vals=None, x_variable=None, y_variable=None,
                  range_x=None, range_y=None, intersections=None, image_widget=None, highlight='', **kwargs):
    colors = None
    image_widget_colorbar =[]
    for i in kwargs:
        p_vals[str(i)] = 10**kwargs[str(i)]
    x_range = range_x
    y_range = range_y

    if image_widget.description == 'Figure_No_Colorbar':
        fig = plt.Figure(figsize=(5, 4), dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.6, 0.7])
    else:
        fig = plt.Figure(figsize=(5, 4), dpi=600, facecolor='w')
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])

    fig = plt.gcf()
    ax = plt.gca()

    if image_widget.description == 'Figure_No_Colorbar':
        colors = ds.draw_2D_slice(ax, p_vals, str(x_variable), str(y_variable),
                                 x_range, y_range,
                         intersections=intersections, colorbar=False)
    else:
        ds.draw_2D_slice(ax, p_vals, str(x_variable), str(y_variable),
                         x_range, y_range,
                         intersections=intersections)


    highlight = str(highlight)
    if highlight != '':
        try:
            case = ds(highlight)
            p_bounds = dict(p_vals)
            p_bounds[x_variable] = x_range
            p_bounds[y_variable] = y_range
            if case.is_valid(p_bounds=p_bounds, strict=False) is True:
                case.draw_2D_slice(ax, p_vals, str(x_variable), str(y_variable),
                                   x_range, y_range, fc='none', ec='k', lw='2.')
        except:
            pass
    ax.plot(log10(p_vals[x_variable]), log10(p_vals[y_variable]), 'k.')
    if colors is not None:
        image_widget_colorbar = colorbar_tabs(colors)
    if image_widget is not None:
        fig = plt.gcf()
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        if image_widget.description != 'Figure_No_Colorbar':
            image_widget.value = data
        else:
            image_ = image_widget.children[0].children[0]
            color_bar = image_widget.children[1]
            image_.value = data
            color_bar.value = image_widget_colorbar

        plt.close()


    return
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_slice_notebook(self, p_vals, x_variable, y_variable, range_x, range_y, slider_ranges, image_container=None, **kwargs):

    # image_container.selected_index = 1
    # image_container.set_title(0, 'Figure')
    # image_container.set_title(1, 'Colorbar')
    # image_container.selected_index = 0

    plot_widget = interactive(make_2D_slice, ds=fixed(self), 
                              p_vals=fixed(p_vals),
                              x_variable=fixed(x_variable),
                              y_variable=fixed(y_variable), 
                              range_x=fixed(range_x),
                              range_y=fixed(range_y),
                              intersections={'Single': [1],
                                             'Triple': [3],
                                             'Single and Triple': [1, 3],
                                             'All': range(1, 100),
                                             },
                              highlight=Text(value=''),
                              image_widget=fixed(image_container),
                              **{i:FloatSlider(min=log10(slider_ranges[i][0]),
                                               max=log10(slider_ranges[i][1]),
                                               step=log10(slider_ranges[i][1]/slider_ranges[i][0])/20,
                                               value=log10(slider_ranges[i][0]*slider_ranges[i][1])/2 if log10(p_vals[i])==0 else log10(p_vals[i]))   #value=log10(p_vals[i]))
                              for i in slider_ranges}
                              )

    make_2D_slice(ds=self, p_vals=p_vals,
                  x_variable=x_variable, y_variable=y_variable,
                  range_x=range_x, range_y=range_y, 
                  image_widget=image_container, #intersections=range(1,100) # intersections=range(1, 100),  for single and triple [1, 3]
                  highlight='', **kwargs)
    return plot_widget
    

