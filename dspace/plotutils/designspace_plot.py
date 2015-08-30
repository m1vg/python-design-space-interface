''' Definition of the abstract model class.


'''

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

from dspace.plotutils.monkey_patching import monkeypatch_method
from dspace.models.designspace import DesignSpace

import dspace.plotutils.case_plot
from dspace.models.designspace import sort_cases

import StringIO
from subprocess import call, Popen, PIPE
from dspace.graphs.designspace_graph import GraphGenerator



lb_plot_colors = {'>,>':(0.75/255, 137/255, 208/255),
                  '0,>':(1., 0.5, 0.),
                  '>,0':(0.5, 0., 1.),
                  '<,>':'g',
                  '>,<':'y',
                  '0,<':(0.95, .95, 0.),
                  '<,0':(0.5, 1., 0.),
                  '<,<':'r',
                  '0,0':'w'}

def key_sort_function(x, y):
    
    if ',' in x and ',' not in y:
        return 1
    if ',' not in x and ',' in y:
        return -1
    if ',' in x and ',' in y:
        x = x.split(',')
        y = y.split(',')
        if len(x) < len(y):
            return -1
        if len(y) < len(x):
            return 1
        for i in xrange(len(x)):
            diff = sort_cases(x, y)
            if diff == 0:
                continue
            return diff
    return sort_cases(x, y)
    return 0
   
 

class SliderCallback(object):
    
    def __init__(self, ds, slider_dictionary, c_axs,  colordict, *args, **kwargs):
        self.ds = ds
        self.sliders = slider_dictionary
        self.c_axs = c_axs
        self.args = args
        self.kwargs = kwargs
        kwargs['color_dict'] = colordict
        kwargs['colorbar'] = False
        
    def __call__(self, val):
        ax = self.args[0]
        pvals = self.args[1]
        sliders = self.sliders
        number_c_axs = len(self.c_axs)
        color_dict = self.kwargs['color_dict']
        for i in sliders:
            pvals[i] = 10**sliders[i].val
        ax.clear()
        for i in self.c_axs:
            i.clear()
        color_dict.update(self.ds.draw_2D_slice(*self.args, **self.kwargs))
        labels = color_dict.keys()
        try:
            labels.sort(cmp=key_sort_function)
        except:
            pass
        labels.reverse()
        num = 0 
        j = 0
        while num < len(labels):
            temp_dict = {i:color_dict[i] for i in labels[num:min(num+20, len(labels))]}
            if j < len(self.c_axs):
                c_ax = self.c_axs[j]
            else:
                c_ax,kw=mt.colorbar.make_axes(ax)
                c_ax.set_aspect(15)
                self.c_axs.append(c_ax)
            self.ds.draw_region_colorbar(c_ax, temp_dict)
            num += 20 
            j += 1
        plt.draw()
             
@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_region_colorbar(self, ax, color_dict, **kwargs):
    
    i = 0
    increment = 1./len(color_dict)
    labels = color_dict.keys()
    try:
        labels.sort(cmp=key_sort_function)
    except:
        pass
    labels.reverse()
    for key in labels:
        ax.fill([0, 1, 1, 0], [i, i, i+increment, i+increment], ec='none', fc=color_dict[key])
        i += increment
    ax.set_yticks([increment*(i+0.5) for i in range(0, len(labels)+1)])
    ax.set_ylim([0, 1])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_ticks_position('right')
    ax.set_yticklabels(labels, **kwargs)
    ax.yaxis.set_ticks_position('none')

@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_function_colorbar(self, ax, zlim, cmap, **kwargs):
    zrange = zlim[1]-zlim[0]
    x=np.linspace(0, 1, 2)
    y=np.linspace(zlim[0], zlim[1], 100)
    a=np.outer(y,np.ones(2))
    ax.pcolor(x, y, a, cmap=cmap, rasterized=True)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ndigits = -int(floor(log10(zrange)))
    attempts = [5, 6, 4, 7, 3, 8, 2]
    ax.set_yticks([round(zlim[0]+i*zrange, ndigits+2) for i in [0., 0.25, 0.5, 0.75, 1.]])
    ax.set_ylim(zlim[0], zlim[1])

@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_log_gain_repertoire(self, ax, x_variable, y_variable, z_variable, 
                                color_dict=lb_plot_colors):
    behavior_set = self.data_2D_log_gain_repertoire(x_variable, 
                                                    y_variable, 
                                                    z_variable)
    for X in behavior_set:
        key = ''
        if X[0] < 0.0:
            key += '<,'
        elif X[0] == 0.0:
            key += '0,'
        elif X[0] > 0.0:
            key += '>,'
        if X[1] < 0.0:
            key += '<'
        elif X[1] == 0.0:
            key += '0'
        elif X[1] > 0.0:
            key += '>'
        C = color_dict[key]
        if X[2] == 0:
            symbol='o'
        else:
            symbol='p'
        ax.plot(X[0], X[1], '.', marker=symbol, mfc=C, mec='k', lw=0.3, ms=5)
    ax.plot([0, 0], [0, -8.5], ls='-', c='gray', lw=0.5)
    ax.plot([0,-8.5], [0, 0], ls='-', c='gray', lw=0.5)
    ax.plot([0, 0], [0, 8.5], ls='-', c='gray', lw=0.5)
    ax.plot([0, 8.5], [0, 0], ls='-', c='gray', lw=0.5)
    ax.set_xticks([-8, -4, 0, 4, 8])
    ax.set_yticks([-8, -4, 0, 4, 8])
    ax.set_xlim([-9, 9])
    ax.set_ylim([-9, 9])

    

@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_routh_index(self, ax, p_vals, x_variable, y_variable, range_x, range_y, color_dict=None,
                           colorbar=True, resolution=100, cmap=mt.cm.Spectral_r):
    
    if color_dict is None:
        color_dict = dict()
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    ssystems = list()
    for case_number in valid_cases:
        case = self(case_number)
        V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                               range_x=[range_x[0]/2, range_x[1]*2], range_y=[range_y[0]/2, range_y[1]*2],
                               log_out=True)
        path=mt.path.Path(V, closed=False) 
        ssystems.append((path, case.ssystem.remove_algebraic_constraints()))
    x = np.linspace(log10(range_x[0]), log10(range_x[1]), resolution)
    y = np.linspace(log10(range_y[0]), log10(range_y[1]), resolution)
    X, Y = np.meshgrid(x, y)
    Z=np.zeros((len(x), len(y)), dtype=np.int)
    values = dict()
    params = VariablePool(p_vals)
    Zj = set()
    for i in xrange(len(y)):
        yi = y[i]
        params[y_variable] = 10**yi
        for j in xrange(len(x)):
            xj = x[j]
            params[x_variable] = 10**xj
            Zj.clear()
            for path, ssystem in ssystems:
                if path.contains_point((xj, yi)):
                    Zj.add(ssystem.routh_index(params))
            nums = [num for num in Zj]
            if len(nums) == 0:
                continue
            key = str(nums[0])
            for index in xrange(1, len(nums)):
                key += ','+str(nums[index])
            try:
                Z[i,j] = values[key]
            except KeyError:
                Z[i,j] = len(values)
                values[key] = len(values)
    colors = dict()
    for i in values:
        colors[values[i]] = cmap(values[i]/(len(values)))
        if i in color_dict:
            colors[values[i]] = color_dict[i]
    fc = [colors[i] for i in xrange(len(values))]
    cf=ax.contourf(X, Y, Z, cmap=None,
                   levels=[-1] + [i for i in xrange(len(values)+1)],
                   colors = fc + ['k'])
    colors = {key:colors[values[key]] for key in values}
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, colors)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    return cf, colors

@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_2D_positive_roots(self, ax, p_vals, x_variable, y_variable, range_x, 
                           range_y, color_dict=None, colorbar=True, 
                           resolution=100, cmap=mt.cm.jet,
                           included_cases=None):
    
    if color_dict is None:
        color_dict = dict()
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    hatched_cases = []
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
    ssystems = list()
    for case_number in valid_cases:
        case = self(case_number)
        if isinstance(case, CyclicalCase) is True:
            V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                                       range_x=[range_x[0]/2, range_x[1]*2],
                                       range_y=[range_y[0]/2, range_y[1]*2],
                                       log_out=True)
            for subcase in V:
                path=mt.path.Path(V[subcase], closed=False) 
                ssystems.append((path, case))
        else:
            V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                                       range_x=[range_x[0]/2, range_x[1]*2],
                                       range_y=[range_y[0]/2, range_y[1]*2],
                                       log_out=True)
            path=mt.path.Path(V, closed=False) 
            ssystems.append((path, case.ssystem.remove_algebraic_constraints()))
    x = np.linspace(log10(range_x[0]), log10(range_x[1]), resolution)
    y = np.linspace(log10(range_y[0]), log10(range_y[1]), resolution)
    X, Y = np.meshgrid(x, y)
    Z=np.zeros((len(x), len(y)), dtype=np.int)
    Z = Z - 1
    values = dict()
    params = VariablePool(p_vals)
    for i in xrange(len(y)):
        yi = y[i]
        params[y_variable] = 10**yi
        for j in xrange(len(x)):
            xj = x[j]
            params[x_variable] = 10**xj
            Zj = []
            for path, system in ssystems:
                if path.contains_point((xj, yi)):
                    roots = system.positive_roots(params)
                    if isinstance(roots, dict) is True:
                        for subcase in roots:
                            Zj.append(roots[subcase])
                    else:
                        Zj.append(roots)
            Zj.sort()
            nums = [num for num in Zj]
            if len(nums) == 0:
                continue
            key = str(nums[0])
            for index in xrange(1, len(nums)):
                key += ','+str(nums[index])
            try:
                Z[i,j] = values[key]
            except KeyError:
                Z[i,j] = len(values)
                values[key] = len(values)
    colors = dict()
    for i in values:
        colors[values[i]] = cmap(values[i]/(len(values)))
        if i in color_dict:
            colors[values[i]] = color_dict[i]
    fc = [colors[i] for i in xrange(len(values))]
    cf=ax.contourf(X, Y, Z, cmap=None,
                   levels=[-2, -1] + [i for i in xrange(len(values)+1)],
                   colors = ['k'] + fc + ['k'])
    colors = {key:colors[values[key]] for key in values}
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, colors)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    return cf, colors
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_slice(self, ax, p_vals, x_variable, y_variable,
                  range_x, range_y, color_dict=None,
                  intersections=[1,2,3,4,5], included_cases=None, 
                  expand_cycles=True,
                  colorbar=True, cmap=mt.cm.gist_rainbow, **kwargs):
    pvals = dspace.VariablePool(names=self.independent_variables)
    if set(pvals.keys()) != set(p_vals.keys()):
        raise ValueError, 'Incomplete parameter set'
    pvals.update(p_vals)
    p_vals = pvals 
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    hatched_cases = []
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_nonstrict = [i for i in valid_cases if i in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
            valid_nonstrict = [i for i in valid_cases if i not in valid_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
        valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
    if len(valid_cases)+len(valid_nonstrict) == 0:
        # fill black
        return
    case_int_list = self.intersecting_cases(intersections, valid_cases+valid_nonstrict, 
                                            p_bounds=p_bounds, strict=False)
    colors = dict()
    if color_dict is None:
        color_dict = dict()
    if 'ec' not in kwargs:
        kwargs['ec']='none'
    ## hatched_cases = self.cycles_to_subcases(hatched_cases)
    for case_num in hatched_cases:
        case = self(case_num)
        case.draw_2D_slice(ax, p_vals, x_variable, y_variable,
                           range_x, range_y, fc='none', 
                           ec=(0.8, 0.8, 0.8, 1.), hatch='/', lw=0.5)
    for case_int in case_int_list:
        key = str(case_int)
        case_nums = key.split(', ')
        for i in xrange(len(case_nums)):
            if expand_cycles is False:
                case_nums[i] = str(case_nums[i]).split('_')[0]
            if case_nums[i] in valid_nonstrict:
                case_nums[i] = str(case_nums[i])+'*'
        key = ', '.join(case_nums)
        if key not in color_dict:
            color_dict[key] = cmap((1.*case_int_list.index(case_int))/len(case_int_list))
        V = case_int.draw_2D_slice(ax, p_vals, x_variable, y_variable,
                                   range_x, range_y, fc=color_dict[key], **kwargs)
        colors[key] = color_dict[key]
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    if colorbar is False:
        return color_dict
    if colorbar is True or colorbar == 'auto':
        labels = colors.keys()
        try:
            labels.sort(cmp=key_sort_function)
        except:
            pass
        labels.reverse()
        num = 0 
        while num < len(labels):
            temp_dict = {i:colors[i] for i in labels[num:min(num+20, len(labels))]}
            c_ax,kw=mt.colorbar.make_axes(ax)
            c_ax.set_aspect(15)
            self.draw_region_colorbar(c_ax, temp_dict)
            num += 20
        plt.sca(ax)
    elif colorbar == 'single':
        labels = colors.keys()
        try:
            labels.sort(cmp=key_sort_function)
        except:
            pass
        labels.reverse()
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, colors)
        plt.sca(ax)
    return color_dict

@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_slice_interactive(self, p_vals, x_variable, y_variable,
                              range_x, range_y, slider_ranges,
                              **kwargs):
    previous = plt.isinteractive()
    plt.ioff()
    number_of_sliders = len(slider_ranges)
    slider_block = 0.03*number_of_sliders
    fig = plt.figure()
    plt.clf()
    ax = plt.axes([0.1, 0.2+slider_block, 0.8, 0.7-slider_block])
    c_axs = list()
    cdict = dict()
    if 'color_dict' in kwargs:
        cdict.update(kwargs['color_dict'])
    j = 0
    sliders = dict()
    for i in slider_ranges:
        slider_ax = plt.axes([0.1, 0.1+j*0.03, 0.8, 0.02])
        slider = Slider(slider_ax, i, 
                        log10(slider_ranges[i][0]), log10(slider_ranges[i][1]), 
                        valinit=log10(p_vals[i]), color='#AAAAAA'
                        )
        j += 1
        sliders[i] = slider
    update = SliderCallback(self, sliders, c_axs, cdict, 
                            ax, p_vals, x_variable, y_variable, range_x, range_y,
                            **kwargs)
    update(1)
    for i in sliders:
        sliders[i].on_changed(update)
    plt.show()
    plt.interactive(previous)
    
 
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_3D_slice(self, ax, p_vals, x_variable, y_variable,z_variable, range_x,
                  range_y, range_z, color_dict=None,
                  intersections=[1,2,3,4,5], included_cases=None, 
                  colorbar=True, cmap=mt.cm.gist_rainbow, **kwargs):
    pvals = dspace.VariablePool(names=self.independent_variables)
    if set(pvals.keys()) != set(p_vals.keys()):
        raise ValueError, 'Incomplete parameter set'
    pvals.update(p_vals)
    p_vals = pvals 
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    p_bounds[z_variable] = range_z
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_nonstrict = [i for i in valid_cases if i in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
            valid_nonstrict = [i for i in valid_cases if i not in valid_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
        valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
    case_int_list = self.intersecting_cases(intersections, 
                                            valid_cases,
                                            p_bounds=p_bounds)
    if color_dict is None:
        color_dict = dict()
    for case_int in case_int_list:
        key = str(case_int)
        if key not in color_dict:
            color_dict[key] = cmap((1.*case_int_list.index(case_int))/len(case_int_list))
        V = case_int.draw_3D_slice(ax, p_vals, x_variable, y_variable, z_variable,
                                   range_x, range_y, range_z, fc=color_dict[key], 
                                   **kwargs)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    ax.set_zlim([log10(min(range_z)), log10(max(range_z))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    if z_variable in self._latex:
        z_variable = '$'+self._latex[z_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    ax.set_zlabel(r'$\log_{10}$(' + z_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, color_dict)
    return color_dict

def calculate_case_2D_function(case, function, p_vals, x_variable, y_variable, 
                        range_x, range_y, resolution, log_linear):
    X, Y, Z, clim, patch = case.draw_2D_ss_function_data(function, 
                                                         p_vals, x_variable, y_variable,
                                                         range_x, range_y, 
                                                         resolution=resolution,
                                                         log_linear=log_linear
                                                         )
    return case, X, Y, Z, clim, patch

@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def _draw_2D_ss_function_parallel(self, ax, all_cases, function, p_vals, x_variable, y_variable, 
                        range_x, range_y, resolution=100, log_linear=False, 
                        zlim=None, included_cases=None, colorbar=True,
                        cmap=mt.cm.jet, surface=False, **kwargs):
    try:
        from joblib import Parallel, delayed, cpu_count
    except:
        return None
    min_lim = 1e20
    max_lim = -1e20
    patches = []
    cores = cpu_count()
    results = Parallel(n_jobs=min(cores, len(all_cases)))(delayed(calculate_case_2D_function)(self(i), 
                                                                      function,
                                                                      p_vals, 
                                                                      x_variable,
                                                                      y_variable,
                                                                      range_x,
                                                                      range_y,
                                                                      resolution,
                                                                      log_linear)
                                                                      for i in all_cases)
    for i in results:
        case, X, Y, Z, clim, path = i
        pc = case.draw_2D_ss_function_from_data(ax, X, Y, Z, clim, path, 
                                                cmap=cmap, 
                                                zlim=zlim, surface=surface,
                                                **kwargs)
        if isinstance(pc, list) is True:
            for apc in pc:
                lims = apc.get_clim()
                min_lim = min(min_lim, lims[0])
                max_lim = max(max_lim, lims[1])
                patches.append(apc)
        else:
            lims = pc.get_clim()
            min_lim = min(min_lim, lims[0])
            max_lim = max(max_lim, lims[1])
            patches.append(pc)
    if zlim is None:
        if min_lim == max_lim:
            delta_z = 1e-3
            min_lim = min_lim-delta_z
            max_lim = max_lim+delta_z
        ndigits = -int(floor(log10(max_lim - min_lim)))
        zlim = [round(min_lim, ndigits), round(max_lim, ndigits)]
    if zlim[0] == zlim[1]:
        delta_z = 1e-3
        zlim = [zlim[0]-delta_z, zlim[1]+delta_z]
    for pc in patches:
        pc.set_clim(zlim)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        self.draw_function_colorbar(c_ax, zlim, cmap)
        c_ax.set_aspect(15./(zlim[1]-zlim[0]))
    return patches
            
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_ss_function(self, ax, function, p_vals, x_variable, y_variable, 
                        range_x, range_y, resolution=100, log_linear=False, 
                        zlim=None, included_cases=None, colorbar=True,
                        cmap=mt.cm.jet, parallel=False, surface=False, **kwargs):                         
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    hatched_cases = []
    valid_nonstrict=[]
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
        valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
    min_lim = 1e20
    max_lim = -1e20
    constant_vars = dict(p_vals)
    constant_vars.pop(x_variable)
    constant_vars.pop(y_variable)
    expr = Expression(function)
    expr = expr.subst(**constant_vars)
    all_cases = list()
    for case in valid_cases+valid_nonstrict:
        if case in valid_nonstrict:
            vertices = self(case).vertices_2D_slice(p_vals, x_variable, y_variable,
                                                    range_x=range_x, range_y=range_y)
            if len(vertices) <= 2:
                continue
            all_cases.append(case)
        else:
            all_cases.append(case)
    if parallel is False:
        patches = None
    else:  
        patches = self._draw_2D_ss_function_parallel(ax,
                                                     all_cases,
                                                     function, 
                                                     p_vals,
                                                     x_variable, 
                                                     y_variable,
                                                     range_x,
                                                     range_y,
                                                     resolution=resolution,
                                                     log_linear=log_linear,
                                                     zlim=zlim,
                                                     included_cases=included_cases,
                                                     colorbar=colorbar,
                                                     cmap=cmap, surface=surface,
                                                     **kwargs)
    if patches is not None:
        return patches  
    patches = list()
    for case in all_cases:
        pc = self(case).draw_2D_ss_function(ax, expr, p_vals, x_variable, y_variable,
                                            range_x, range_y, resolution=resolution,
                                            log_linear=log_linear, cmap=cmap, 
                                            surface=surface, **kwargs)
        if isinstance(pc, list) is True:
            for apc in pc:
                lims = apc.get_clim()
                min_lim = min(min_lim, lims[0])
                max_lim = max(max_lim, lims[1])
                patches.append(apc)
        else:
            lims = pc.get_clim()
            min_lim = min(min_lim, lims[0])
            max_lim = max(max_lim, lims[1])
            patches.append(pc)
    if zlim is None:
        if min_lim == max_lim:
            delta_z = 1e-3
            min_lim = min_lim-delta_z
            max_lim = max_lim+delta_z
        ndigits = -int(floor(log10(max_lim - min_lim)))
        zlim = [round(min_lim, ndigits), round(max_lim, ndigits)]
    if zlim[0] == zlim[1]:
        delta_z = 1e-3
        zlim = [zlim[0]-delta_z, zlim[1]+delta_z]
    for pc in patches:
        pc.set_clim(zlim)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        self.draw_function_colorbar(c_ax, zlim, cmap)
        c_ax.set_aspect(15./(zlim[1]-zlim[0]))
    return patches

@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_dominant_eigenvalues(self, ax, p_vals, x_variable, y_variable, 
                                 range_x, range_y, resolution=100, component='real', 
                                 zlim=None, included_cases=None, colorbar=True,
                                 cmap=mt.cm.jet, parallel=False, cmp=None, **kwargs):                         
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    hatched_cases = []
    valid_nonstrict=[]
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
        valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
    min_lim = 1e20
    max_lim = -1e20
    all_cases = list()
    for case in valid_cases+valid_nonstrict:
        if case in valid_nonstrict:
            vertices = self(case).vertices_2D_slice(p_vals, x_variable, y_variable,
                                                    range_x=range_x, range_y=range_y)
            if len(vertices) <= 2:
                continue
            all_cases.append(case)
        else:
            all_cases.append(case)
    if parallel is False:
        patches = None
    else:  
        patches = None
        ## self._draw_2D_ss_function_parallel(ax,
        ##                                              all_cases,
        ##                                              function, 
        ##                                              p_vals,
        ##                                              x_variable, 
        ##                                              y_variable,
        ##                                              range_x,
        ##                                              range_y,
        ##                                              resolution=resolution,
        ##                                              log_linear=log_linear,
        ##                                              zlim=zlim,
        ##                                              included_cases=included_cases,
        ##                                              colorbar=colorbar,
        ##                                              cmap=cmap, **kwargs)
    if patches is not None:
        return patches  
    patches = list()
    for case in all_cases:
        pc = self(case).draw_2D_dominant_eigenvalue(ax, p_vals, 
                                                    x_variable, y_variable,
                                                    range_x, range_y,
                                                    resolution=resolution,
                                                    component=component,
                                                    cmap=cmap, cmp=cmp,
                                                    **kwargs)
        if isinstance(pc, list) is True:
            for apc in pc:
                lims = apc.get_clim()
                min_lim = min(min_lim, lims[0])
                max_lim = max(max_lim, lims[1])
                patches.append(apc)
        else:
            lims = pc.get_clim()
            min_lim = min(min_lim, lims[0])
            max_lim = max(max_lim, lims[1])
            patches.append(pc)
    if zlim is None:
        if min_lim == max_lim:
            delta_z = 1e-3
            min_lim = min_lim-delta_z
            max_lim = max_lim+delta_z
        ndigits = -int(floor(log10(max_lim - min_lim)))
        zlim = [round(min_lim, ndigits), round(max_lim, ndigits)]
    if zlim[0] == zlim[1]:
        delta_z = 1e-3
        zlim = [zlim[0]-delta_z, zlim[1]+delta_z]
    for pc in patches:
        pc.set_clim(zlim)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    if x_variable in self._latex:
        x_variable = '$'+self._latex[x_variable]+'$'
    if y_variable in self._latex:
        y_variable = '$'+self._latex[y_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        self.draw_function_colorbar(c_ax, zlim, cmap)
        c_ax.set_aspect(15./(zlim[1]-zlim[0]))
    return patches
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_1D_slice(self, ax, p_vals, slice_variable, range_slice, color_dict=None,
                  intersections=[1,2,3,4,5], colorbar=True, cmap=mt.cm.gist_rainbow,
                  included_cases = None, **kwargs):
    
    p_bounds = dict(p_vals)
    p_bounds[slice_variable] = range_slice
    ## valid_cases = self.valid_cases(p_bounds=p_bounds)
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_nonstrict = [i for i in valid_cases if i in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
            valid_nonstrict = [i for i in valid_cases if i not in valid_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
        valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
    if len(valid_cases)+len(valid_nonstrict) == 0:
        # fill black
        return
    case_int_list = self.intersecting_cases(intersections, valid_cases+valid_nonstrict, 
                                            p_bounds=p_bounds, strict=False)
    if color_dict is None:
        color_dict = dict()
        
    for case_int in case_int_list:
        key = str(case_int)
        if key not in color_dict:
            color_dict[key] = cmap((1.*case_int_list.index(case_int))/len(case_int_list))
        V = case_int.draw_1D_slice(ax, p_vals, slice_variable, range_slice, fc=color_dict[key],
                                   **kwargs)
    ax.set_xlim([log10(min(range_slice)), log10(max(range_slice))])
    if slice_variable in self._latex:
        slice_variable = '$'+self._latex[slice_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + slice_variable + ')')
    ax.yaxis.set_visible(False)
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, color_dict)
    return color_dict
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_1D_ss_function(self, ax, function, p_vals, 
                        slice_variable, range_slice, 
                        resolution=100, colors=None, included_cases=None, ylim = None, **kwargs):
    p_bounds = dict(p_vals)
    p_bounds[slice_variable] = range_slice
    if included_cases is not None:
        included_cases = [i.case_number for i in self(included_cases)]
        if self.number_of_cases < 1e5:
            valid_cases = self.valid_cases(p_bounds=p_bounds)
            valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
            hatched_cases = [i for i in valid_cases if i not in included_cases]
            valid_nonstrict = [i for i in valid_cases if i in included_cases]
            valid_cases = [i for i in valid_cases if i in included_cases]
            valid_nonstrict = [i for i in valid_cases if i not in valid_cases]
        else:
            valid_cases = [i for i in included_cases if self(i).is_valid(p_bounds=p_bounds)]
    else:
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        valid_nonstrict = self.valid_cases(p_bounds=p_bounds, strict=False)
        valid_nonstrict = [i for i in valid_nonstrict if i not in valid_cases]
    if len(valid_cases)+len(valid_nonstrict) == 0:
        # fill black
        return
    valid_cases = valid_cases + valid_nonstrict#self.valid_cases(p_bounds=p_bounds)
    lines = list()
    ylim_t = None
    if 'color' not in kwargs:
        kwargs['color'] = 'k'
    for case in valid_cases:
        if colors is not None:
            if case in colors:
                kwargs['color'] = colors[case]
            else:
                kwargs.pop('color')
        pt = self(case).draw_1D_ss_function(ax, 
                                            function,
                                            p_vals,
                                            slice_variable, 
                                            range_slice,
                                            resolution=resolution,
                                            **kwargs)
        if pt is None:
            continue
        lines.append(pt)
        ydata = pt[0].get_ydata()
        miny = min(ydata)
        maxy = max(ydata)
        if ylim is None:
            if ylim_t is None:
                ylim_t = [min(ydata), max(ydata)]
            else:
                ylim_t = [min((ylim_t[0], miny)), max((ylim_t[1], maxy))]
    if ylim is None:
        ylim = ylim_t
    ax.set_ylim([ylim[0]-(ylim[1]-ylim[0])*0.1, ylim[1]+(ylim[1]-ylim[0])*0.1])
    ax.set_xlim(np.log10(range_slice))
    if slice_variable in self._latex:
        slice_variable = '$'+self._latex[slice_variable]+'$'
    ax.set_xlabel(r'$\log_{10}$(' + slice_variable + ')')
    if isinstance(function, Expression):
        expr = function
    else:
        expr = Expression(function)
    ax.set_ylabel('$'+expr.latex(self._latex)+'$')    
    return lines
   
@monkeypatch_method(dspace.models.designspace.DesignSpace) 
def draw_1D_positive_roots(self, ax, function, p_vals, slice_variable, 
                           range_slice, resolution=100, ylim=None,
                           line_dict=None, **kwargs):
    lines = self.line_1D_positive_roots(function, p_vals, slice_variable,
                                        range_slice, resolution=int(resolution))
    line_styles = ['-', '--', '..']
    colors = ['k', 'r', 'y']
    unique_R = {i[2] for i in lines}
    if line_dict == None:
        line_dict = {0:{'ls':'-','c':'k','lw':2.},
                     1:{'ls':':','c':'r','lw':2.},
                     2:{'ls':'--','c':'y','lw':2.}
                     }
        j = 0
        for i in unique_R:
            if i in line_dict:
                continue
            line_dict[i] = {'c':mt.cm.Spectral_r(j/len(unique_R)),
                            'ls':'--'}
    for i in lines:
        x = i[0]
        y = i[1]
        r = i[2]
        ax.plot(x, y, **line_dict[r])
    if slice_variable in self._latex:
        slice_variable = '$'+self._latex[slice_variable]+'$'
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(r'$\log_{10}$(' + slice_variable + ')')
    if isinstance(function, Expression):
        expr = function
    else:
        expr = Expression(function)
    ax.set_ylabel('$'+expr.latex(self._latex)+'$')
    return lines

@monkeypatch_method(dspace.models.designspace.DesignSpace) 
def draw_network_graph(self, ax, graph_type='dot', p_vals=None, 
                       included_variables=[], size=[3.33, 2],
                       resolution=100, cmap=mt.cm.jet, colorbar=True,
                       show_regulation=True,
                       **kwargs):
    g = GraphGenerator(self)
    dot_data = g.graph(graph_type=graph_type,
                       p_vals=p_vals,
                       cmap=cmap,
                       included_variables=included_variables,
                       show_regulation=show_regulation)
    dot_string = dot_data['description']
    cmd = Popen(['dot', '-Tpng', '-Gdpi='+str(resolution)],
                stdin=PIPE, stdout=PIPE, stderr=PIPE)
    data, err = cmd.communicate(input=dot_string)
    if len(err) > 0:
        raise OSError, err
    if 'limits' in dot_data:
        if colorbar is True:
            c_ax,kw=mt.colorbar.make_axes(ax)
            zlim = dot_data['limits']
            self.draw_function_colorbar(c_ax, zlim, cmap)
            c_ax.set_aspect(15./(zlim[1]-zlim[0]))
    f=StringIO.StringIO(data)
    cax=ax.imshow(mt.image.imread(f))
    ax.xaxis.visible=False
    ax.yaxis.visible=False
    ax.set_xticks([])
    ax.set_yticks([])
    plt.box(False)
    f.close()
