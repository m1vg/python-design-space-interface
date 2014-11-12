''' Definition of the abstract model class.


'''

import numpy as np
import matplotlib as mt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as a3
from math import log10

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
import dspace.models.case

from dspace.plotutils.monkey_patching import monkeypatch_method
from dspace.models.case import Case, CaseIntersection


def generate_plot_lattice_bounds(case, p_vals, x_variable, y_variable, range_x, range_y, resolution):
    
    V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                               range_x=range_x, range_y=range_y,
                               log_out=True)
    path=mt.path.Path(V)    
    bbox=path.get_extents()
    fraction_x = (bbox.max[0]-bbox.min[0])/(log10(range_x[1])-log10(range_x[0]))
    fraction_y = (bbox.max[1]-bbox.min[1])/(log10(range_y[1])-log10(range_y[0]))
    resolution_x = int(resolution*fraction_x)+5
    resolution_y = int(resolution*fraction_y)+5
    x=np.linspace(bbox.min[0], bbox.max[0], resolution_x+5)
    y=np.linspace(bbox.min[1], bbox.max[1], resolution_y+5)
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

@monkeypatch_method(dspace.models.case.Case)
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
   
    
def interpolated_data(case, function, p_vals, x_variable, y_variable, points, x, y, path):

    f_val = list()
    params = VariablePool(names=case.independent_variables)
    for key in params:
        params[key] = p_vals[key]
    clim = None
    for (x_value,y_value) in points:
        params[x_variable] = 10**x_value
        params[y_variable] = 10**y_value
        value = case.ssystem.steady_state_function(function, params)        
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
    return X,Y,Z,clim

def sampled_data(case, function, p_vals, x_variable, y_variable, points, x, y, path):

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
            Z[i,j] = case.ssystem.steady_state_function(function, params)
            if path.contains_point((xj, yi)) == 1:
                if clim is None:
                    clim = [Z[i,j], Z[i,j]]
                else:
                    clim[0] = min(clim[0], Z[i,j])
                    clim[1] = max(clim[1], Z[i,j])
    return X,Y,Z,clim
         
@monkeypatch_method(dspace.models.case.Case)
def draw_2D_ss_function_data(self, function, p_vals, x_variable, y_variable,
                        range_x, range_y, 
                        resolution=100, log_linear=False): 
    points, x, y, path = generate_plot_lattice_bounds(self, p_vals,
                                                      x_variable, y_variable,
                                                      range_x, range_y, resolution)
    if log_linear is True:
        X,Y,Z,clim = interpolated_data(self, function, p_vals,
                                  x_variable, y_variable,
                                  points, x, y, path)
    else:
        X,Y,Z,clim = sampled_data(self, function, p_vals,
                             x_variable, y_variable,
                             points, x, y, path)
    ## patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    return (X, Y, Z, clim, path)
    pc = ax.pcolor(X, Y, Z, rasterized=True, **kwargs)
    ax.add_patch(patch)
    pc.set_clip_path(patch)
    if zlim is None:
        zlim=clim
    if zlim is not None:    
        pc.set_clim(zlim)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return pc
        
    
@monkeypatch_method(dspace.models.case.Case)
def draw_2D_ss_function_from_data(self, ax, X, Y, Z, clim, path, zlim=None, **kwargs):
    patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    pc = ax.pcolor(X, Y, Z, rasterized=True, **kwargs)
    ax.add_patch(patch)
    pc.set_clip_path(patch)
    if zlim is None:
        zlim=clim
    if zlim is not None:    
        pc.set_clim(zlim)
    return pc
         
@monkeypatch_method(dspace.models.case.Case)
def draw_2D_ss_function(self, ax, function, p_vals, x_variable, y_variable,
                        range_x, range_y, 
                        resolution=100, log_linear=False, zlim=None, **kwargs):
    
    X, Y, Z, clim, path = self.draw_2D_ss_function_data(function, 
                                                         p_vals,
                                                         x_variable,
                                                         y_variable,
                                                         range_x,
                                                         range_y,
                                                         resolution=resolution,
                                                         log_linear=log_linear
                                                         )
    ## patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    pc = self.draw_2D_ss_function_from_data(ax, X, Y, Z, clim, path, zlim=zlim, **kwargs)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return pc

@monkeypatch_method([dspace.models.case.Case, dspace.models.case.CaseIntersection])   
def draw_2D_slice(self, ax, p_vals, x_variable, y_variable, range_x, range_y,
                  show_equations=False, rotation=30, fontsize=8, **kwargs):
    
    if show_equations is False:
        V = self.vertices_2D_slice(p_vals, x_variable, y_variable,
                                   range_x=range_x, range_y=range_y,
                                   log_out=True)
    else:
        vertices = self.vertices_2D_slice(p_vals, x_variable, y_variable,
                                          range_x=range_x, range_y=range_y,
                                          log_out=False, vtype='both')
        V = [(log10(i[0][0]), log10(i[0][1])) for i in vertices]
    if len(V) == 0:
        return
    V = zip(*V)
    ax.fill(V[0], V[1], **kwargs)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    if show_equations is True:
        for i in xrange(len(vertices)):
            x = log10(vertices[i][0][0])
            y = log10(vertices[i][0][1])
            s = '\n'.join(['$'+j.latex(self._latex)+'$' for j in vertices[i][1]])
            ax.plot(x, y, 'k.', mfc='none', mec='k', ms=1.)
            ax.text(x, y, s, fontsize=fontsize, rotation=rotation, 
                    horizontalalignment='center', verticalalignment='center')

@monkeypatch_method([dspace.models.case.Case, dspace.models.case.CaseIntersection])   
def draw_3D_slice(self, ax, p_vals, x_variable, y_variable, z_variable, range_x, range_y, range_z,
                  **kwargs):
    
    faces = self.faces_3D_slice(p_vals, x_variable, y_variable, z_variable,
                             range_x=range_x, range_y=range_y, range_z=range_z,
                             log_out=True)
    tri = a3.art3d.Line3DCollection(faces)
    if 'fc' in kwargs:
        tri.set_facecolor(kwargs['fc'])
    elif 'facecolor' in kwargs:
        tri.set_facecolor(kwargs['facecolor'])
    else:
        tri.set_facecolor('none')
    tri.set_alpha(0.5)
    if 'ec' in kwargs:
        tri.set_edgecolor(kwargs['ec'])
    elif 'edgecolor' in kwargs:
        tri.set_edgecolor(kwargs['edgecolor'])
    else:
        tri.set_edgecolor('k')
    ax.add_collection3d(tri)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    ax.set_zlim(np.log10(range_z))

    
@monkeypatch_method(dspace.models.case.Case)
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
    ax.set_xlim(np.log10(range_slice))
    print self
    return pt

@monkeypatch_method([dspace.models.case.Case, dspace.models.case.CaseIntersection])
def draw_1D_slice(self, ax, p_vals, slice_variable, range_slice, **kwargs):
    
    V = self.vertices_1D_slice(p_vals, slice_variable, range_slice=range_slice, log_out=True)
    V = zip(*V)[0]
    pt = ax.fill([V[0], V[1], V[1], V[0]], [0, 0, 1, 1], **kwargs)
    #pt = ax.plot(V, [0, 0], **kwargs)
    return pt
