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


def make_2D_slice(ds=None, p_vals=None, x_variable=None, y_variable=None,
                  range_x=None, range_y=None, **kwargs):
    for i in kwargs:
#         display(str(i))
        p_vals[str(i)] = 10**kwargs[str(i)]
    x_range = range_x
    y_range = range_y
    ax = plt.gca()
    ds.draw_2D_slice(ax, p_vals, str(x_variable), str(y_variable),
                             x_range, y_range)
    ax.plot(log10(pvals[x_variable]), log10(pvals[y_variable]), 'k.')

@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_slice_notebook(self, p_vals, x_variable, y_variable,
                           range_x, range_y, slider_ranges,
                           **kwargs):
    plot_widget = interact(make_2D_slice, ds=fixed(self), 
                           p_vals=fixed(p_vals), 
                           x_variable=fixed(x_variable),
                           y_variable=fixed(y_variable), 
                           range_x=fixed(range_x),
                           range_y=fixed(range_y),
                           **{i:widgets.FloatSliderWidget(min=log10(slider_ranges[i][0]),
                                                          max=log10(slider_ranges[i][1]),
                                                          step=log10(slider_ranges[i][1]/slider_ranges[i][0])/20,
                                                          value=log10(p_vals[i]))
                              for i in slider_ranges})
    