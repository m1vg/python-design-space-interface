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
    resolution_x = int(resolution*fraction_x)+2
    resolution_y = int(resolution*fraction_y)+2
    x=np.linspace(bbox.min[0], bbox.max[0], resolution_x+2)
    y=np.linspace(bbox.min[1], bbox.max[1], resolution_y+2)
    points = [bbox.min, (bbox.min[0], bbox.max[1]), bbox.max, (bbox.max[0], bbox.min[1])]
    return (points,x,y,path)

def generate_plot_lattice_bounds_new(case, p_vals, x_variable, y_variable, range_x, range_y, resolution):
    
    V = case.vertices_2D_slice(p_vals, x_variable, y_variable, 
                               range_x=range_x, range_y=range_y,
                               log_out=True)
    path=mt.path.Path(V)    
    bbox=path.get_extents()
    ## fraction_x = (bbox.max[0]-bbox.min[0])/(log10(range_x[1])-log10(range_x[0]))
    ## fraction_y = (bbox.max[1]-bbox.min[1])/(log10(range_y[1])-log10(range_y[0]))
    ## resolution_x = int(resolution*fraction_x)+2
    ## resolution_y = int(resolution*fraction_y)+2
    
    x=np.linspace(log10(range_x[0]), log10(range_x[1]), resolution)
    y=np.linspace(log10(range_y[0]), log10(range_y[1]), resolution)
    i= 0
    j = 0
    x_indices = []
    y_indices = []
    for xi in x[:-1]:
        if xi > bbox.min[0]:
            break
        i += 1
    x_indices.append(max(0, i-1))
    for xi in x[i:-1]:
        if xi > bbox.max[0]:
            break
        i += 1    
    x_indices.append(min(len(y), i+1))
    for yi in y[:-1]:
        if yi > bbox.min[1]:
            break
        j += 1
    y_indices.append(max(0, j-1))
    for yi in y[j:-1]:
        if yi > bbox.max[1]:
            break
        j += 1   
    y_indices.append(min(len(y), j+1)) 
         
    return (x_indices,y_indices,path)

                
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
   
    
def interpolated_data(case, function, p_vals, x_variable, y_variable, points, x, y, path, alt_function=None):

    f_val = list()
    alt_f_val = list()
    params = VariablePool(names=case.independent_variables)
    for key in params:
        params[key] = p_vals[key]
    clim = None
    for (x_value,y_value) in points:
        params[x_variable] = 10**x_value
        params[y_variable] = 10**y_value
        value = case.ssystem.steady_state_function(function, params)        
        f_val.append(value)
        if alt_function is not None:
            alt_value = case.ssystem.steady_state_function(alt_function, params)
            alt_f_val.append(alt_value)
    V = zip(*points)
    X,Y = np.meshgrid(x, y)
    Z = mt.mlab.griddata(V[0], V[1], f_val, X, Y, interp='linear')
    if alt_function is not None:
        Z_alt = mt.mlab.griddata(V[0], V[1], alt_f_val, X, Y, interp='linear')
    for (i, yi) in enumerate(y):
        for (j, xj) in enumerate(x):
            if path.contains_point((xj, yi)) == 1:
                if clim is None:
                    clim = [Z[i,j], Z[i,j]]
                else:
                    clim[0] = min(clim[0], Z[i,j])
                    clim[1] = max(clim[1], Z[i,j])
    if alt_function is not None:
        Z = [Z, Z_alt]
    return X,Y,Z,clim

def sampled_data(case, function, p_vals, x_variable, y_variable, points, x, y, path, alt_function=None):

    f_val = list()
    alt_f_val = list()
    params = VariablePool(p_vals)
    V = zip(*points)    
    X,Y = np.meshgrid(x, y)
    Z = mt.mlab.griddata(V[0], V[1], np.repeat(1, len(V[0])), X, Y)
    if alt_function is not None:
        Z_alt = mt.mlab.griddata(V[0], V[1], np.repeat(1, len(V[0])), X, Y)
    clim = None
    for (i, yi) in enumerate(y):
        params[y_variable] = 10**y[i]
        for (j, xj) in enumerate(x):
            params[x_variable] = 10**x[j]
            Z[i,j] = case.ssystem.steady_state_function(function, params)
            if alt_function is not None:
                Z_alt[i,j] = case.ssystem.steady_state_function(alt_function, params)
            if path.contains_point((xj, yi)) == 1:
                if clim is None:
                    clim = [Z[i,j], Z[i,j]]
                else:
                    clim[0] = min(clim[0], Z[i,j])
                    clim[1] = max(clim[1], Z[i,j])
    if alt_function is not None:
        Z = [Z, Z_alt]
    return X,Y,Z,clim

def interpolated_data_new(case, function, p_vals, x_variable, y_variable,
                          range_x, range_y, x_indices, y_indices, path,
                          resolution, alt_function=None):
    f_val = list()
    alt_f_val = list()
    params = VariablePool(names=case.independent_variables)
    for key in params:
        params[key] = p_vals[key]
    clim = None
    points = path.vertices
    for (x_value,y_value) in points:
        params[x_variable] = 10**x_value
        params[y_variable] = 10**y_value
        value = case.ssystem.steady_state_function(function, params)        
        f_val.append(value)
        if alt_function is not None:
            alt_value = case.ssystem.steady_state_function(alt_function, params)
            alt_f_val.append(alt_value)
    V = zip(*points)
    delta_x = (range_x[1]-range_x[0])/resolution
    delta_y = (range_y[1]-range_y[0])/resolution
    x = np.linspace(range_x[0] + delta_x*x_indices[0], range_x[0] + delta_x * x_indices[1], 1+x_indices[1] - x_indices[0])
    y = np.linspace(range_y[0] + delta_y*y_indices[0], range_y[0] + delta_y * y_indices[1], 1+y_indices[1] - y_indices[0])
    x = np.round(x, 10)
    y = np.round(y, 10)
    X,Y = np.meshgrid(x, y)
    Z = mt.mlab.griddata(V[0], V[1], f_val, X, Y, interp='linear')
    if alt_function is not None:
        Z_alt = mt.mlab.griddata(V[0], V[1], alt_f_val, X, Y, interp='linear')
    for (i, yi) in enumerate(y):
        for (j, xj) in enumerate(x):
            if path.contains_point((xj, yi)) == False:
                Z[i,j] = np.nan
                continue
            if clim is None:
                clim = [Z[i,j], Z[i,j]]
            else:
                clim[0] = min(clim[0], Z[i,j])
                clim[1] = max(clim[1], Z[i,j])
    Z = np.ma.array(Z, mask=np.isnan(Z))
    interpolated_path = path.interpolated((len(path.vertices)-1)*resolution)
    for (xj, yi) in interpolated_path.vertices:
        j = np.argmin(abs(x-xj))
        i = np.argmin(abs(y-yi))
        params[x_variable] = 10**x[j]
        params[y_variable] = 10**y[i]
        Z[i,j] = case.ssystem.steady_state_function(function, params)
        if alt_function is not None:
            Z_alt[i,j] = case.ssystem.steady_state_function(alt_function, params)
        if clim is None:
            clim = [Z[i,j], Z[i,j]]
        else:
            clim[0] = min(clim[0], Z[i,j])
            clim[1] = max(clim[1], Z[i,j])
    if alt_function is not None:
        Z = [Z, Z_alt]
    return X,Y,Z,clim

def sample_data_new(case, function, p_vals, x_variable, y_variable,
                     range_x, range_y, x_indices, y_indices, path, resolution, alt_function=None):
    f_val = list()
    alt_f_val = list()
    
    params = VariablePool(p_vals)  
    delta_x = (range_x[1]-range_x[0])/resolution
    delta_y = (range_y[1]-range_y[0])/resolution
    x = np.linspace(range_x[0] + delta_x*x_indices[0], 
                    range_x[0] + delta_x * x_indices[1], 
                    1+x_indices[1] - x_indices[0])
    y = np.linspace(range_y[0] + delta_y*y_indices[0], 
                    range_y[0] + delta_y * y_indices[1], 
                    1+y_indices[1] - y_indices[0])
    X,Y = np.meshgrid(x, y)
    Z = np.zeros((len(y), len(x)))
    if alt_function is not None:
        Z = np.zeros((len(x), len(y))) #mt.mlab.griddata(x, y, np.repeat(0, len(x)), X, Y)
    clim = None
    for (i, yi) in enumerate(y):
        params[y_variable] = 10**yi
        for (j, xj) in enumerate(x):
            if path.contains_point((xj, yi)) == False:
                Z[i,j] = np.nan
                continue
            params[x_variable] = 10**xj
            Z[i,j] = case.ssystem.steady_state_function(function, params)
            if alt_function is not None:
                Z_alt[i,j] = case.ssystem.steady_state_function(alt_function, params)
            if clim is None:
                clim = [Z[i,j], Z[i,j]]
            else:
                clim[0] = min(clim[0], Z[i,j])
                clim[1] = max(clim[1], Z[i,j])
    Z = np.ma.array(Z, mask=np.isnan(Z))
    interpolated_path = path.interpolated((len(path.vertices)-1)*resolution)
    for (xj, yi) in interpolated_path.vertices:
        j = np.argmin(abs(x-xj))
        i = np.argmin(abs(y-yi))
        params[x_variable] = 10**x[j]
        params[y_variable] = 10**y[i]
        ## print xj, yi
        Z[i,j] = case.ssystem.steady_state_function(function, params)
        if alt_function is not None:
            Z_alt[i,j] = case.ssystem.steady_state_function(alt_function, params)
        if clim is None:
            clim = [Z[i,j], Z[i,j]]
        else:
            clim[0] = min(clim[0], Z[i,j])
            clim[1] = max(clim[1], Z[i,j])
    if alt_function is not None:
        Z = [Z, Z_alt]
    return X,Y,Z,clim
    

@monkeypatch_method(dspace.models.case.Case)
def draw_2D_dominant_eigenvalue_data(self, p_vals, x_variable, y_variable,
                                     range_x, range_y, cmp=None,
                                     resolution=100, component='real'): 
    params = VariablePool(p_vals)
    clim = None  
    x_indices, y_indices, path = generate_plot_lattice_bounds_new(self, p_vals,
                                                  x_variable, y_variable,
                                                  range_x, range_y, resolution)
    range_x = np.log10(range_x)
    range_y = np.log10(range_y)
    delta_x = (range_x[1]-range_x[0])/resolution
    delta_y = (range_y[1]-range_y[0])/resolution
    x = np.linspace(range_x[0] + delta_x*x_indices[0], 
                    range_x[0] + delta_x * x_indices[1], 
                    1+x_indices[1] - x_indices[0])
    y = np.linspace(range_y[0] + delta_y*y_indices[0], 
                    range_y[0] + delta_y * y_indices[1], 
                    1+y_indices[1] - y_indices[0])
    x = np.round(x, 10)
    y = np.round(y, 10)
    X,Y = np.meshgrid(x, y)
    Z = np.zeros((len(y), len(x)))
    ssys = self.ssystem.remove_algebraic_constraints()
    for (i, yi) in enumerate(y):
        params[y_variable] = 10**yi
        for (j, xj) in enumerate(x):
            if path.contains_point((xj, yi)) == False:
                Z[i,j] = np.nan
                continue
            params[x_variable] = 10**xj
            eigen_values = ssys.eigenvalues(params)
            if cmp is None:
                if component == 'real':
                    Z[i,j] = max(eigen_values.real)
                else:
                    Z[i,j] = max(eigen_values.imag)
            else:
                Z[i,j] = cmp(eigen_values)
            if clim is None:
                clim = [Z[i,j], Z[i,j]]
            else:
                clim[0] = min(clim[0], Z[i,j])
                clim[1] = max(clim[1], Z[i,j])
    Z = np.ma.array(Z, mask=np.isnan(Z))
    interpolated_path = path.interpolated((len(path.vertices)-1)*resolution)
    for (xj, yi) in interpolated_path.vertices:
        j = np.argmin(abs(x-xj))
        i = np.argmin(abs(y-yi))
        params[x_variable] = 10**x[j]
        params[y_variable] = 10**y[i]
        eigen_values = ssys.eigenvalues(params)
        if cmp is None:
            if component == 'real':
                Z[i,j] = max(eigen_values.real)
            else:
                Z[i,j] = max(eigen_values.imag)
        else:
            Z[i,j] = cmp(eigen_values)
        if clim is None:
            clim = [Z[i,j], Z[i,j]]
        else:
            clim[0] = min(clim[0], Z[i,j])
            clim[1] = max(clim[1], Z[i,j])
    return (X, Y, Z, clim, path)
        
@monkeypatch_method(dspace.models.case.Case)
def draw_2D_dominant_eigenvalue(self, ax, p_vals, x_variable, y_variable,
                                range_x, range_y, resolution=100,
                                zlim=None, component='real', cmp=None,
                                **kwargs):
    
    X, Y, Z, clim, path = self.draw_2D_dominant_eigenvalue_data(p_vals,
                                                                x_variable,
                                                                y_variable,
                                                                range_x,
                                                                range_y,
                                                                resolution=resolution,
                                                                component=component,
                                                                cmp=cmp,
                                                                )
    if 'cmap' in kwargs:
        cmap = kwargs.pop('cmap') 
        cmap.set_bad((0, 0, 0, 0))
    else:
        cmap = mt.cm.jet
        cmap.set_bad((0., 0., 0., 0.))
    pc = self.draw_2D_ss_function_from_data(ax, X, Y, Z, clim, path,
                                            zlim=zlim, cmap=cmap, **kwargs)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return pc

@monkeypatch_method(dspace.models.case.Case)
def draw_2D_ss_function_data(self, function, p_vals, x_variable, y_variable,
                             range_x, range_y, 
                             resolution=100, log_linear=False): 
    ## points, x, y, path = generate_plot_lattice_bounds(self, p_vals,
    ##                                                   x_variable, y_variable,
    ##                                                   range_x, range_y, resolution)
    x, y, path = generate_plot_lattice_bounds_new(self, p_vals,
                                                          x_variable, y_variable,
                                                          range_x, range_y, resolution)
    if log_linear is True:
        X,Y,Z,clim = interpolated_data_new(self, function, p_vals,
                                           x_variable, y_variable,
                                           [log10(i) for i in range_x], 
                                           [log10(i) for i in range_y],
                                           x, y, path, resolution)
    else:
        X,Y,Z,clim = sample_data_new(self, function, p_vals,
                                     x_variable, y_variable,
                                     [log10(i) for i in range_x], 
                                     [log10(i) for i in range_y],
                                     x, y, path, resolution)
    ## patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    return (X, Y, Z, clim, path)        
    
@monkeypatch_method(dspace.models.case.Case)
def draw_2D_ss_function_from_data(self, ax, X, Y, Z, clim, path, zlim=None, surface=False, **kwargs):
    patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=0.)
    if surface is True:
        pc = ax.plot_surface(X, Y, Z, edgecolor='none')
    else:
        pc = ax.pcolormesh(X, Y, Z, rasterized=True, **kwargs)
    ## ax.add_patch(patch)
    ## pc.set_clip_path(patch)
    if zlim is None:
        zlim=clim
    if zlim is not None:    
        pc.set_clim(zlim)
    return pc

    
@monkeypatch_method(dspace.models.case.Case)
def draw_2D_ss_function(self, ax, function, p_vals, x_variable, y_variable,
                        range_x, range_y, resolution=100, log_linear=False, zlim=None, **kwargs):
    
    X, Y, Z, clim, path = self.draw_2D_ss_function_data(function, 
                                                         p_vals,
                                                         x_variable,
                                                         y_variable,
                                                         range_x,
                                                         range_y,
                                                         resolution=resolution,
                                                         log_linear=log_linear
                                                         )
    if 'cmap' in kwargs:
        cmap = kwargs.pop('cmap') 
        cmap.set_bad((0, 0, 0, 0))
    else:
        cmap = mt.cm.jet
        cmap.set_bad((0., 0., 0., 0.))
    ## patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    pc = self.draw_2D_ss_function_from_data(ax, X, Y, Z, clim, path, zlim=zlim, cmap=cmap, **kwargs)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return pc


@monkeypatch_method(dspace.models.case.Case)
def draw_2D_phase_portrait_data(self, p_vals, x_variable, y_variable,
                                range_x, range_y, x_dynamic, y_dynamic, 
                                resolution=100, log_linear=False): 
    points, x, y, path = generate_plot_lattice_bounds(self, p_vals,
                                                      x_variable, y_variable,
                                                      range_x, range_y, resolution)
    if x_dynamic == True:
        function = 'log($e_'+x_variable+'_p) - log($e_'+x_variable+'_n)'
        ## function = x_variable+'*(log($e_'+x_variable+'_p) - log($e_'+x_variable+'_n))'
        ## function = '$e_'+x_variable+'_p - $e_'+x_variable+'_n'
    else:
        function = '0'
    if y_dynamic == True:
        alt_function = 'log($e_'+y_variable+'_p) - log($e_'+y_variable+'_n)'
        ## alt_function = y_variable+'*(log($e_'+y_variable+'_p) - log($e_'+y_variable+'_n))'
        ## alt_function = '$e_'+y_variable+'_p - $e_'+y_variable+'_n'
    else:
        alt_function = '0'
    if log_linear is True:
        X,Y,Z,clim = interpolated_data(self, function, p_vals,
                                  x_variable, y_variable,
                                  points, x, y, path, alt_function=alt_function)
    else:
        X,Y,Z,clim = sampled_data(self, function, p_vals,
                             x_variable, y_variable,
                             points, x, y, path, alt_function=alt_function)
    ## patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    Zx = Z[0]
    Zy = Z[1]
    return (X, Y, Zx, Zy, path)
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
def draw_2D_phase_portrait_from_data(self, ax, X, Y, Zx, Zy, path, log_linear=False, **kwargs):
    patch = mt.patches.PathPatch(path, fc='none', ec='none', lw=.5)
    q = ax.quiver(X, Y, Zx, Zy, angles='xy',
                  #scale_units='xy', 
                  #scale=3, 
                  headwidth=2., 
                  headlength=2, 
                  headaxislength=2.25,
                  pivot='mid', **kwargs)
    ax.add_patch(patch)
    q.set_clip_path(patch)
    return q

@monkeypatch_method(dspace.models.case.Case)
def draw_2D_phase_portrait(self, ax, Xd_initial, p_vals, x_variable, y_variable,
                        range_x, range_y, 
                        resolution=100, log_linear=False, zlim=None, show_designspaces=False, color_dict=None, **kwargs):
    x_dynamic = False
    y_dynamic = False
    if x_variable in self.dependent_variables:
        x_dynamic = True
    if y_variable in self.dependent_variables:
        y_dynamic=True
    if x_dynamic is False and y_dynamic is False:
        return None
    es = self.eigen_spaces()
    p = VariablePool(names=es.independent_variables)
    for i in p:
        if i in Xd_initial:
            p[i] = Xd_initial[i]
        elif i in p_vals:
            p[i] = p_vals[i]
    bounds = dict(p)
    bounds[x_variable] = range_x
    bounds[y_variable] = range_y
    valid = es.valid_cases(p_bounds=bounds, strict=False)
    Q = list()
    if show_designspaces is True:
        es.draw_2D_slice(ax, p, x_variable, y_variable,
                         range_x, range_y, 
                         colorbar=False, color_dict=color_dict,
                         alpha=0.2)
    for i in valid:
        space = es(i)
        vertices = space.vertices_2D_slice(p, x_variable, y_variable,
                                           range_x=range_x, range_y=range_y, log_out=True)
        if len(vertices) <= 2:
            continue
        X, Y, Zx, Zy, path = space.draw_2D_phase_portrait_data(p, 
                                                               x_variable,
                                                               y_variable,
                                                               range_x,
                                                               range_y,
                                                               x_dynamic,
                                                               y_dynamic,
                                                               resolution=resolution,
                                                               log_linear=log_linear)
        q = self.draw_2D_phase_portrait_from_data(ax, X, Y, Zx, Zy, path, **kwargs)
        Q.append(q)
    ax.set_xlim(np.log10(range_x))
    ax.set_ylim(np.log10(range_y))
    return Q

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
    if len(V) == 0:
        return None
    V = zip(*V)[0]
    if len(V) == 1:
        return None
    X = np.linspace(V[0], V[1], resolution)
    f_val = list()
    for x in X:
        params[slice_variable] = 10**x
        f_val.append(self.ssystem.steady_state_function(function, params))
    pt = ax.plot(X, f_val, **kwargs)
    ax.set_xlim(np.log10(range_slice))
    return pt

@monkeypatch_method([dspace.models.case.Case, dspace.models.case.CaseIntersection])
def draw_1D_slice(self, ax, p_vals, slice_variable, range_slice, **kwargs):
    
    V = self.vertices_1D_slice(p_vals, slice_variable, range_slice=range_slice, log_out=True)
    if len(V) == 0:
        return None
    V = zip(*V)[0]
    if len(V) == 1:
        return None
    pt = ax.fill([V[0], V[1], V[1], V[0]], [0, 0, 1, 1], **kwargs)
    #pt = ax.plot(V, [0, 0], **kwargs)
    return pt
