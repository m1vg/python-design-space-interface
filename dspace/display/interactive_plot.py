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


def make_2D_slice(ds=None, p_vals=None, x_variable=None, y_variable=None,
                  range_x=None, range_y=None, intersections=None, image_widget=None, highlight='', **kwargs):
    for i in kwargs:
        p_vals[str(i)] = 10**kwargs[str(i)]
    x_range = range_x
    y_range = range_y
    fig = plt.Figure(figsize=[6, 4], dpi=600, facecolor='w')
    fig=plt.gcf()
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax = plt.gca()
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
    if image_widget is not None:
        fig = plt.gcf()
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        image_widget.value = data
        plt.close()
    return
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_slice_notebook(self, p_vals, x_variable, y_variable,
                           range_x, range_y, slider_ranges,
                           image_container=None, **kwargs):
    plot_widget = interactive(make_2D_slice, ds=fixed(self), 
                              p_vals=fixed(p_vals),
                              x_variable=fixed(x_variable),
                              y_variable=fixed(y_variable), 
                              range_x=fixed(range_x),
                              range_y=fixed(range_y),
                              intersections={'Single':[1],
                                             'Single and Triple':[1,3],
                                             'Triple':[3],
                                             'All':range(1, 100)},
                              highlight=Text(value=''),
                              image_widget=fixed(image_container),
                              **{i:FloatSlider(min=log10(slider_ranges[i][0]),
                                                             max=log10(slider_ranges[i][1]),
                                                             step=log10(slider_ranges[i][1]/slider_ranges[i][0])/20,
                                                             value=log10(p_vals[i]))
                                                             for i in slider_ranges
                                                             }

                                 )
    make_2D_slice(ds=self, p_vals=p_vals,
                  x_variable=x_variable, y_variable=y_variable,
                  range_x=range_x, range_y=range_y, 
                  intersections=range(1, 100), image_widget=image_container, 
                  highlight='')
    return plot_widget
    