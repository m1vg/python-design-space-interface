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

from IPython.html.widgets import interact, interactive, fixed
from IPython.html import widgets
from IPython.display import clear_output, display, HTML

import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg  


def make_2D_slice(ds=None, p_vals=None, x_variable=None, y_variable=None,
                  range_x=None, range_y=None, intersections=None, image_widget=None, **kwargs):
    for i in kwargs:
        p_vals[str(i)] = 10**kwargs[str(i)]
    x_range = range_x
    y_range = range_y
    fig = plt.Figure(figsize=[6, 6])
    fig=plt.gcf()
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
    ax = plt.gca()
    ds.draw_2D_slice(ax, p_vals, str(x_variable), str(y_variable),
                             x_range, y_range,
                     intersections=intersections)
    ax.plot(log10(p_vals[x_variable]), log10(p_vals[y_variable]), 'k.')
    if image_widget is not None:
        fig = plt.gcf()
        canvas = FigureCanvasAgg(fig) 
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        data = buf.getvalue()
        image_widget.value = data
        fig.clf()
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
                                             'All':range(1, 100)},
                              image_widget=fixed(image_container),
                              **{i:widgets.FloatSliderWidget(min=log10(slider_ranges[i][0]),
                                                             max=log10(slider_ranges[i][1]),
                                                             step=log10(slider_ranges[i][1]/slider_ranges[i][0])/20,
                                                          value=log10(p_vals[i]))
                                                          for i in slider_ranges
                                                          }
                              )
    return plot_widget
    