from __future__ import division

import numpy as np
import matplotlib as mt
import matplotlib.pyplot as plt
from math import log10

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
import dspace.models.case

from dspace.plotutils.monkey_patching import monkeypatch_method
from dspace.models.case import Case, CaseIntersection

from dspace.plotutils import case_plot


def generate_plot_lattice_bounds(case, p_vals, x_variable, y_variable, range_x, range_y, resolution):
    
    V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                               range_x=range_x, range_y=range_y,
                               log_out=True)
    path=mt.path.Path(V)    
    bbox=path.get_extents()
    fraction_x = (bbox.max[0]-bbox.min[0])/(log10(range_x[1])-log10(range_x[0]))
    fraction_y = (bbox.max[1]-bbox.min[1])/(log10(range_y[1])-log10(range_y[0]))
    resolution_x = int(resolution*fraction_x)+2
    resolution_y = int(resolution*fraction_y)+2
    x=np.linspace(bbox.min[0], bbox.max[0], resolution_x)
    y=np.linspace(bbox.min[1], bbox.max[1], resolution_y)
    points = [bbox.min, (bbox.min[0], bbox.max[1]), bbox.max, (bbox.max[0], bbox.min[1])]
    return (points,x,y,path)
                
def routh_data(case, p_vals, x_variable, y_variable, points, x, y, path):

    f_val = list()
    params = VariablePool(p_vals)
    V = zip(*points)    
    X,Y = np.meshgrid(x, y)
    Z = mt.mlab.griddata(V[0], V[1], np.repeat(1, len(V[0])), X, Y)
    levels = list()
    ssystem = case.ssystem.remove_algebraic_constraints()
    mask=np.ones((len(y), len(x)))
    for (i, yi) in enumerate(y):
        params[y_variable] = 10**y[i]
        for (j, xj) in enumerate(x):
            params[x_variable] = 10**x[j]
            Z[i,j] = ssystem.routh_index(params)
            if path.contains_point((xj, yi)) == 1:
                if Z[i,j] not in levels:
                    levels.append(Z[i,j])
                mask[i,j] = 0
    Z = np.ma.masked_array(Z, mask=mask)
    return X,Y,Z,levels

@monkeypatch_method(dspace.models.cyclicalcase.CyclicalCase)
def draw_2D_routh(self, ax, p_vals, x_variable, y_variable, range_x, range_y, 
                  resolution=100, log_linear=False, levels=None, **kwargs):
    
    points, x, y, path = generate_plot_lattice_bounds(self, p_vals,
                                                      x_variable, y_variable,
                                                      range_x, range_y, resolution)    
    X,Y,Z,lvls = routh_data(self, p_vals, x_variable, y_variable,
                            points, x, y, path)
    patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    if len(lvls) is 0:
        lvls = [0]
    if levels is None:
        levels = lvls
        draw_levels = list(levels)
        draw_levels.append(max(draw_levels)+1)
    else:
        draw_levels = set(levels)
    cf = ax.contourf(X, Y, Z, rasterized=True, levels=draw_levels, **kwargs)
    ax.add_patch(patch)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return cf,levels
   
    
def interpolated_data(cyclicalcase, subcase, function, p_vals, x_variable, y_variable, points, x, y, path):

    f_val = list()
    params = VariablePool(names=cyclicalcase.independent_variables)
    for key in params:
        params[key] = p_vals[key]
    clim = None
    for (x_value,y_value) in points:
        params[x_variable] = 10**x_value
        params[y_variable] = 10**y_value
        f_vals = cyclicalcase.steady_state_function(function, params)
        value = f_vals[str(subcase)]
        f_val.append(value)
    V = zip(*points)
    X,Y = np.meshgrid(x, y)
    Z = mt.mlab.griddata(V[0], V[1], f_val, X, Y, interp='linear')
    for (i, yi) in enumerate(y):
        for (j, xj) in enumerate(x):
            if path.contains_point((xj, yi)) == 1:
                if clim is None:
                    clim = [Z[i,j], Z[i,j]]
                else:
                    clim[0] = min(clim[0], Z[i,j])
                    clim[1] = max(clim[1], Z[i,j])
##     return X,Y,Z,clim

def sampled_data(cyclicalcase, subcase, function, p_vals, x_variable, y_variable, points, x, y, path):

    f_val = list()
    params = VariablePool(p_vals)
    V = zip(*points)    
    X,Y = np.meshgrid(x, y)
    Z = mt.mlab.griddata(V[0], V[1], np.repeat(1, len(V[0])), X, Y)
    clim = None
    for (i, yi) in enumerate(y):
        params[y_variable] = 10**y[i]
        for (j, xj) in enumerate(x):
            params[x_variable] = 10**x[j]
            ## print cyclicalcase.steady_state_function(function, params)
            try:
                Z[i,j] = cyclicalcase.steady_state_function(function, params)[str(subcase)]
                if path.contains_point((xj, yi)) == 1:
                    if clim is None:
                        clim = [Z[i,j], Z[i,j]]
                    else:
                        clim[0] = min(clim[0], Z[i,j])
                        clim[1] = max(clim[1], Z[i,j])
            except:
                pass
    return X,Y,Z,clim
                     
@monkeypatch_method(dspace.models.cyclicalcase.CyclicalCase)
def draw_2D_ss_function(self, ax, function, p_vals, x_variable, y_variable,
                        range_x, range_y, 
                        resolution=100, log_linear=False, zlim=None, **kwargs):
    
    p_bounds = dict(p_vals)
    p_bounds[x_variable] = range_x
    p_bounds[y_variable] = range_y
    subcases = self.valid_subcases(p_bounds=p_bounds)
    pcs = list()
    clim=None

    for subcase in subcases:
        case = self(subcase)
        points, x, y, path = case_plot.generate_plot_lattice_bounds(case, p_vals,
                                                                    x_variable, y_variable,
                                                                    range_x, range_y, resolution)
        if log_linear is True:
            X,Y,Z,clim_sub = interpolated_data(self, subcase, function, p_vals,
                                           x_variable, y_variable,
                                           points, x, y, path)
        else:
            X,Y,Z,clim_sub = sampled_data(self, subcase, function, p_vals,
                                          x_variable, y_variable,
                                          points, x, y, path)
        patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
        if clim == None:
            clim = clim_sub
        else:
            clim = [min([clim[0], clim_sub[0]]),
                    max([clim[1], clim_sub[1]])]
        pc = ax.pcolor(X, Y, Z, rasterized=True, **kwargs)
        ax.add_patch(patch)
        pc.set_clip_path(patch)
        pcs.append(pc)
    if zlim is None:
        zlim=clim
    if zlim is not None:
        for pc in pcs:
            pc.set_clim(zlim)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return pcs

@monkeypatch_method(dspace.models.cyclicalcase.CyclicalCase)   
def draw_2D_slice(self, ax, p_vals, x_variable, y_variable, range_x, range_y,
                  **kwargs):
    
    subcases = self.vertices_2D_slice(p_vals, x_variable, y_variable,
                                      range_x=range_x, range_y=range_y,
                                      log_out=True)
    for case in subcases:
        V = zip(*subcases[case])
        ax.fill(V[0], V[1], **kwargs)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    
@monkeypatch_method(dspace.models.cyclicalcase.CyclicalCase)
def draw_1D_ss_function(self, ax, function, p_vals, slice_variable, range_slice,
                        resolution=100, **kwargs):
    
    params = VariablePool(p_vals)
    V = self.vertices_1D_slice(params, slice_variable, range_slice=range_slice, log_out=True)
    V = zip(*V)[0]
    X = np.linspace(V[0], V[1], resolution)
    f_val = list()
    for x in X:
        params[slice_variable] = 10**x
        f_val.append(self.ssystem.steady_state_function(function, params))
    pt = ax.plot(X, f_val, **kwargs)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return pt

@monkeypatch_method(dspace.models.cyclicalcase.CyclicalCase)
def draw_1D_slice(self, ax, p_vals, slice_variable, range_slice, **kwargs):
    
    V = self.vertices_1D_slice(p_vals, slice_variable, range_slice=range_slice, log_out=True)
    V = zip(*V)[0]
    pt = ax.fill([V[0], V[1], V[1], V[0]], [0, 0, 1, 1], **kwargs)
    #pt = ax.plot(V, [0, 0], **kwargs)
    return pt