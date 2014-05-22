''' Definition of the abstract model class.


'''

from __future__ import division
import numpy as np
import matplotlib as mt
import matplotlib.pyplot as plt
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

def key_sort_function(x, y):
    
    x = x.split(',')
    y = y.split(',')
    if len(x) < len(y):
        return -1
    if len(y) < len(x):
        return 1
    for i in xrange(len(x)):
        if str(x[i]) < str(y[i]):
            return -1
        if str(y[i]) < str(x[i]):
            return 1
    return 0
        
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
    ax.set_yticklabels(labels, **kwargs)
    ax.set_ylim([0, 1])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_ticks_position('right')
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
    ax.set_yticks([round(zlim[0]+i*zrange, 2) for i in [0., 0.25, 0.5, 0.75, 1.]])
    ax.set_ylim(zlim[0], zlim[1])

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
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    if included_cases is not None:
        included_cases = [str(i) for i in included_cases]
        valid_cases = [i for i in valid_cases if str(i) in included_cases]
    ssystems = list()
    for case_number in valid_cases:
        case = self(case_number)
        if isinstance(case, CyclicalCase) is True:
            V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                                       range_x=[range_x[0]/2, range_x[1]*2], range_y=[range_y[0]/2, range_y[1]*2],
                                       log_out=True)
            for subcase in V:
                path=mt.path.Path(V[subcase], closed=False) 
                ssystems.append((path, case))
        else:
            V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                                       range_x=[range_x[0]/2, range_x[1]*2], range_y=[range_y[0]/2, range_y[1]*2],
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
    Zj = set()
    for i in xrange(len(y)):
        yi = y[i]
        params[y_variable] = 10**yi
        for j in xrange(len(x)):
            xj = x[j]
            params[x_variable] = 10**xj
            Zj.clear()
            for path, system in ssystems:
                if path.contains_point((xj, yi)):
                    roots = system.positive_roots(params)
                    if isinstance(roots, dict) is True:
                        for subcase in roots:
                            Zj.add(roots[subcase])
                    else:
                        Zj.add(roots)
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
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    return cf, colors
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_slice(self, ax, p_vals, x_variable, y_variable,
                   range_x, range_y, color_dict=None,
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
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    valid_cases = self.cycles_to_subcases(valid_cases)
    hatched_cases = []
    if included_cases is not None:
        included_cases = [str(i) for i in included_cases]
        hatched_cases = [i for i in valid_cases if str(i) not in included_cases]
        valid_cases = [i for i in valid_cases if str(i) in included_cases]
    case_int_list = self.intersecting_cases(intersections, valid_cases, p_bounds=p_bounds)
    colors = dict()
    if color_dict is None:
        color_dict = dict()
    if 'ec' not in kwargs:
        kwargs['ec']='none'
    hatched_cases = self.cycles_to_subcases(hatched_cases)
    for case_num in hatched_cases:
        case = self(case_num)
        case.draw_2D_slice(ax, p_vals, x_variable, y_variable,
                           range_x, range_y, fc='none', 
                           ec=(0.8, 0.8, 0.8, 1.), hatch='/', lw=0.5)
    for case_int in case_int_list:
        key = str(case_int)
        if key not in color_dict:
            color_dict[key] = cmap((1.*case_int_list.index(case_int))/len(case_int_list))
        V = case_int.draw_2D_slice(ax, p_vals, x_variable, y_variable,
                                   range_x, range_y, fc=color_dict[key], **kwargs)
        colors[key] = color_dict[key]
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, colors)
    return color_dict

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
        included_cases = [str(i) for i in included_cases]
        valid_cases = [i for i in valid_cases if str(i) in included_cases]
    case_int_list = self.intersecting_cases(intersections, 
                                            valid_cases,
                                            p_bounds=p_bounds)
    print intersections
    print case_int_list
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
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    ax.set_zlabel(r'$\log_{10}$(' + z_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, color_dict)
    return color_dict
        
@monkeypatch_method(dspace.models.designspace.DesignSpace)   
def draw_2D_ss_function(self, ax, function, p_vals, x_variable, y_variable, 
                        range_x, range_y, resolution=100, log_linear=False, 
                        zlim=None, included_cases=None, colorbar=True,
                        cmap=mt.cm.jet, **kwargs):
    
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    hatched_cases = []
    if included_cases is not None:
        included_cases = [str(i) for i in included_cases]
        valid_cases = [i for i in valid_cases if str(i) in included_cases]
        hatched_cases = [i for i in valid_cases if str(i) in included_cases]
    min_lim = 1e20
    max_lim = -1e20
    constant_vars = dict(p_vals)
    constant_vars.pop(x_variable)
    constant_vars.pop(y_variable)
    expr = Expression(function)
    expr = expr.subst(**constant_vars)
    patches = list()
    for case in valid_cases:
        pc = self(case).draw_2D_ss_function(ax, expr, p_vals, x_variable, y_variable,
                                            range_x, range_y, resolution=resolution,
                                            log_linear=log_linear, cmap=cmap, **kwargs)
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
        zlim = [floor(min_lim), ceil(max_lim)]
    if zlim[0] == zlim[1]:
        delta_z = zlim[0]*0.1
        zlim = [zlim[0]-delta_z, zlim[1]+delta_z]
    for pc in patches:
        pc.set_clim(zlim)
    ax.set_xlim([log10(min(range_x)), log10(max(range_x))])
    ax.set_ylim([log10(min(range_y)), log10(max(range_y))])
    ax.set_xlabel(r'$\log_{10}$(' + x_variable + ')')
    ax.set_ylabel(r'$\log_{10}$(' + y_variable + ')')
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        self.draw_function_colorbar(c_ax, zlim, cmap)
        c_ax.set_aspect(15./(zlim[1]-zlim[0]))
    return patches

@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_1D_slice(self, ax, p_vals, slice_variable, range_slice, color_dict=None,
                  intersections=[1,2,3,4,5], colorbar=True, cmap=mt.cm.gist_rainbow, **kwargs):
    
    p_bounds = dict(p_vals)
    p_bounds[slice_variable] = range_slice
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    case_int_list = self.intersecting_cases(intersections, valid_cases, p_bounds=p_bounds)
    if color_dict is None:
        color_dict = dict()
        
    for case_int in case_int_list:
        key = str(case_int)
        if key not in color_dict:
            color_dict[key] = cmap((1.*case_int_list.index(case_int))/len(case_int_list))
        V = case_int.draw_1D_slice(ax, p_vals, slice_variable, range_slice, fc=color_dict[key],
                                   **kwargs)
    ax.set_xlim([log10(min(range_slice)), log10(max(range_slice))])
    ax.set_xlabel(r'$\log_{10}$(' + slice_variable + ')')
    ax.yaxis.set_visible(False)
    if colorbar is True:
        c_ax,kw=mt.colorbar.make_axes(ax)
        c_ax.set_aspect(15)
        self.draw_region_colorbar(c_ax, color_dict)
    return color_dict
    
@monkeypatch_method(dspace.models.designspace.DesignSpace)
def draw_1D_ss_function(self, ax, function, p_vals, slice_variable, range_slice,
                        resolution=100, **kwargs):
    p_bounds = dict(p_vals)
    p_bounds[slice_variable] = range_slice
    valid_cases = self.valid_cases(p_bounds=p_bounds)
    lines = list()
    ylim = None
    for case in valid_cases:
        pt = self(case).draw_1D_ss_function(ax, function, p_vals,
                                            slice_variable, range_slice,
                                            resolution=resolution, **kwargs)
        lines.append(pt)
        ydata = pt[0].get_ydata()
        miny = min(ydata)
        maxy = max(ydata)
        if ylim is None:
            ylim = [min(ydata), max(ydata)]
        else:
            ylim = [min((ylim[0], miny)), max((ylim[1], maxy))]
    ax.set_ylim([ylim[0]-(ylim[1]-ylim[0])*0.1, ylim[1]+(ylim[1]-ylim[0])*0.1])
    ax.set_xlim(np.log10(range_slice))
    return lines
   
@monkeypatch_method(dspace.models.designspace.DesignSpace) 
def draw_1D_positive_roots(self, ax, function, p_vals, slice_variable, 
                           range_slice, resolution=100,
                           line_dict=None, **kwargs):
    lines = self.line_1D_positive_roots(function, p_vals, slice_variable,
                                        range_slice, resolution=resolution)
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
    return lines
