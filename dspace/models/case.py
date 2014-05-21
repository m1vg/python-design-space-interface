''' Definition of the abstract model class.


'''

from __future__ import division
from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.models.ssystem import SSystem
from math import *


class Case(Model):
    
    def __init__(self, model, swigwrapper, name=None, constraints=None, **kwargs):
        ''' Init method for the model base class.
        
        The model object is initialized with data to construct
        a dictionary of key value pairs, where the keys represent
        names of the model's dependent variables, and the values
        are the equations in the system.
        '''
        if name == None:
            name = 'Unnamed'
        super(Case, self).__init__(model._equations,
                                   name=name)
        setattr(self, '_ssystem', None)
        setattr(self, '_independent_variables', None)
        self.set_swigwrapper(swigwrapper)
        if constraints is not None:
            DSCaseAddConstraints(self._swigwrapper, constraints, len(constraints))
            
    def __del__(self):
        ''' 
        '''
        if self._swigwrapper is not None:
            self._ssystem.set_swigwrapper(None)
            DSCaseFree(self._swigwrapper)
        
    def __str__(self):
        case_info = self.name.split('Case')[1]
        case_info = case_info.split('Subcase')
        if len(case_info) == 1:
            return case_info[0].strip()
        else:
            return '_'.join([i.strip(' :') for i in case_info])
    def set_swigwrapper(self, case_swigwrapper):
        self._swigwrapper = case_swigwrapper
        Xd = VariablePool()
        Xd.set_swigwrapper(DSVariablePoolCopy(DSSSystemXd(DSCaseSSystem(case_swigwrapper))))
        for i in VariablePool():
            if i not in self.dependent_variables:
                raise NameError, 'Dependent Variables are inconsistent'
        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSSSystemXi(DSCaseSSystem(case_swigwrapper))))
        self._independent_variables = Xi
        eqs = list()
        eqs_expr = DSSSystemEquations(DSCaseSSystem(case_swigwrapper))
        for i in xrange(0, len(super(Case, self).equations)):
            eqs.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        self._ssystem = SSystem(self._equations,
                                name=self.name,
                                swigwrapper=DSCaseSSystem(case_swigwrapper))

    @property
    def equations(self):
        eqs = DSCaseEquations(self._swigwrapper)
        equations = list()
        for i in xrange(0, DSCaseNumberOfEquations(self._swigwrapper)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equations.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        return equations
    
    @property
    def ssystem(self):
        return self._ssystem
    
    @property
    def independent_variables(self):
        return self._independent_variables.keys()
    
    def _valid_parameter_set_bounded(self, p_bounds):
        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        for i in lower:
            lower[i] = 1e-20
            upper[i] = 1e20
            if i not in p_bounds:
                continue
            k = p_bounds[i]
            try:
                lower[i] = min(k)
                upper[i] = max(k)
            except:
                lower[i] = k
                upper[i] = k
        variablepool = DSCaseValidParameterSetAtSlice(self._swigwrapper, 
                                                      lower._swigwrapper,
                                                      upper._swigwrapper)
        pvals = VariablePool()
        pvals.set_swigwrapper(variablepool)
        return pvals
    
    def valid_parameter_set(self, p_bounds=None):
        variablepool = DSCaseValidParameterSet(self._swigwrapper)
        if p_bounds is not None:
            pvals = self._valid_parameter_set_bounded(p_bounds)
            return pvals
        pvals = VariablePool()
        pvals.set_swigwrapper(variablepool)
        return pvals
    
    def valid_interior_parameter_set(self, p_bounds=None, distance=50):
        pvals = self.valid_parameter_set(p_bounds=p_bounds)
        for j in xrange(0, 2):
            for i in pvals:
                if p_bounds is not None:
                    if i not in p_bounds:
                        continue
                    if isinstance(p_bounds[i], list) is False:
                        continue
                    slice_range = self.vertices_1D_slice(pvals, i, 
                                                         range_slice=[max([pvals[i]*(2*distance)**-1, min(p_bounds[i])]), 
                                                                      min([pvals[i]*2*distance, max(p_bounds[i])])],
                                                         log_out=True
                                                         )
                    slice_range = [x[0] for x in slice_range]
                    pvals[i] = 10**(slice_range[0] + (slice_range[1] - slice_range[0])/2)
        return pvals
        
    @property
    def case_number(self):
        return DSCaseNumber(self._swigwrapper)

    @property
    def signature(self):
        return DSCaseSignatureToString(self._swigwrapper)
    
    @property
    def case_signature(self):
        return DSCaseSignatureToString(self._swigwrapper)
    @property
    def conditions(self):
        conditions = list()
        eqs_expr = DSCaseConditions(self._swigwrapper)
        for i in xrange(0, DSCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return conditions
        
    @property
    def conditions_log(self):
        conditions = list()
        eqs_expr = DSCaseLogarithmicConditions(self._swigwrapper)
        for i in xrange(0, DSCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return conditions
    
    @property
    def boundaries(self):
        boundaries = list()
        eqs_expr = DSCaseBoundaries(self._swigwrapper)
        for i in xrange(0, DSCaseNumberOfConditions(self._swigwrapper)):
            boundaries.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return boundaries
    
    @property
    def boundaries_log(self):
        boundaries = list()
        eqs_expr = DSCaseLogarithmicBoundaries(self._swigwrapper)
        for i in xrange(0, DSCaseNumberOfBoundaries(self._swigwrapper)):
            boundaries.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return boundaries
    
    @property
    def is_cyclical(self):
        return False
    
    def _is_valid_slice(self, p_bounds):

        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        for i in lower:
            lower[i] = 1E-20
            upper[i] = 1E20
        for (key,value) in p_bounds.iteritems():
            try:
                min_value,max_value = value
            except TypeError:
                min_value = value
                max_value = value
            if min_value > max_value:
                raise ValueError, 'parameter slice bounds are inverted: min is larger than max'
            lower[key] = min_value
            upper[key] = max_value
        return DSCaseIsValidAtSlice(self._swigwrapper,
                                    lower._swigwrapper,
                                    upper._swigwrapper)
    
    def steady_state(self, parameter_values):
        return self.ssystem.steady_state(parameter_values)
        
    def steady_state_flux(self, parameter_values):
        return self.ssystem.steady_state_flux(parameter_values)
        
    def steady_state_function(self, function, parameter_values):
        return self.ssystem.steady_state_function(function, parameter_values)
    
    def is_consistent(self, p_bounds=None):
        
        if p_bounds is not None:
            raise NotImplementedError, 'Needs to be implemented'
            ## return self._is_valid_slice(p_bounds)
            #do something
        return DSCaseConditionsAreValid(self._swigwrapper)
 
    def is_valid(self, p_bounds=None):
        
        if p_bounds is not None:
            return self._is_valid_slice(p_bounds)
            #do something
        return DSCaseIsValid(self._swigwrapper)
        
    def _is_valid_point_in_statespace(self, v_bounds, p_bounds):

        Xd_t = VariablePool()
        independent = VariablePool(names=self.independent_variables)
        ssys = self.ssystem
        for key in independent:
            value = float(p_bounds[key])
            independent[key] = value
        for key in self.dependent_variables:
            if key in ssys.auxiliary_variables:
                continue
            value = float(v_bounds[key])
            Xd_t[key] = value      
        dependent = VariablePool(names=self.dependent_variables)
        aux = ssys.value_for_auxiliary_variables(Xd_t, independent)
        dependent.update(v_bounds)
        dependent.update(aux)
        return DSCaseIsValidInStateSpaceAtPoint(self._swigwrapper,
                                                dependent._swigwrapper,
                                                independent._swigwrapper)
        
    def is_valid_in_state_space(self, v_bounds=None, p_bounds=None):
        
        if p_bounds is not None:
            return self._is_valid_point_in_statespace(v_bounds, p_bounds)
            #do something
        return DSCaseIsValidInStateSpace(self._swigwrapper)
    
    def bounding_box(self, p_bounds=None, log_out=False):
        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        boundkeys = lower.keys()
        for key in lower:
            lower[key] = 1e-20
            upper[key] = 1e20
        if p_bounds is not None:
            for key,value in p_bounds.iteritems():
                if key not in lower:
                    continue;
                try:
                    min_value,max_value = value
                except TypeError:
                    min_value = value
                    max_value = value
                    boundkeys.remove(key)
                if min_value > max_value:
                    raise ValueError, 'Min cannot be larger than max'
                lower[key] = min_value
                upper[key] = max_value
        box = {}
        for key in boundkeys:
            box[key] = DSCaseBoundingRangeForVariableWithConstraints(self._swigwrapper,
                                                                     key,
                                                                     lower._swigwrapper,
                                                                     upper._swigwrapper)
            if log_out is False:
                box[key] = [10**i[0] for i in box[key]]
            else:
                box[key] = [i[0] for i in box[key]]
        return box
            
    def measure_tolerance(self, pvals, log_out=False):
        tolerances={}
        for key in self.independent_variables:
            tolerances[key] = self.vertices_1D_slice(pvals, key, log_out=log_out)
            if log_out is False:
                tolerances[key] = (tolerances[key][0][0]/pvals[key], tolerances[key][1][0]/pvals[key])
            if log_out is True:
                tolerances[key] = (tolerances[key][0][0]-log10(pvals[key]), tolerances[key][1][0]-log10(pvals[key]))
        return tolerances
              
    def vertices_1D_slice(self, p_vals, slice_variable, range_slice=None, log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        if range_slice is None:
            lower[slice_variable] = 1E-20
            upper[slice_variable] = 1E20
        else:
            lower[slice_variable] = min(range_slice)
            upper[slice_variable] = max(range_slice)
        log_vertices=DSCaseVerticesFor1DSlice(self._swigwrapper, 
                                              lower._swigwrapper,
                                              upper._swigwrapper,
                                              slice_variable)
        vertices = list()
        for vertex in log_vertices:
            vertices.append([10**coordinate for coordinate in vertex])
        if log_out is True:
            return log_vertices
        return vertices
    
    def line_1D_positive_roots(self, function, p_vals, slice_variable, range_slice,
                           resolution=100, **kwargs):
            
            params = VariablePool(p_vals)
            V = self.vertices_1D_slice(params, slice_variable, range_slice=range_slice, log_out=True)
            V = zip(*V)[0]
            step =(V[1]-V[0])/resolution
            X = [V[0]+step*i for i in range(resolution+1)]
            f_val = list()
            roots = list()
            ssys = self.ssystem.remove_algebraic_constraints()
            for x in X:
                params[slice_variable] = 10**x
                f_val.append(self.ssystem.steady_state_function(function, params))
                roots.append(ssys.positive_roots(params))
            return (X, f_val, roots)
            
    def vertices_2D_slice(self, p_vals, x_variable, y_variable, range_x=None, range_y=None,
                       log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        if range_x is None:
            lower[x_variable] = 1E-20
            upper[x_variable] = 1E20
        else:
            lower[x_variable] = min(range_x)
            upper[x_variable] = max(range_x)
        if range_y is None:
            lower[y_variable] = 1E-20
            upper[y_variable] = 1E20
        else:
            lower[y_variable] = min(range_y)
            upper[y_variable] = max(range_y)
        log_vertices=DSCaseVerticesFor2DSlice(self._swigwrapper, 
                                              lower._swigwrapper,
                                              upper._swigwrapper,
                                              x_variable,
                                              y_variable)
        vertices = list()
        for vertex in log_vertices:
            vertices.append([10**coordinate for coordinate in vertex])
        if log_out is True:
            return log_vertices
        return vertices

    def vertices_3D_slice(self, p_vals, x_variable, y_variable, z_variable, 
                          range_x=None, range_y=None, range_z=None,
                          log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        if range_x is None:
            lower[x_variable] = 1E-20
            upper[x_variable] = 1E20
        else:
            lower[x_variable] = min(range_x)
            upper[x_variable] = max(range_x)
        if range_y is None:
            lower[y_variable] = 1E-20
            upper[y_variable] = 1E20
        else:
            lower[y_variable] = min(range_y)
            upper[y_variable] = max(range_y)
        if range_y is None:
            lower[z_variable] = 1E-20
            upper[z_variable] = 1E20
        else:
            lower[z_variable] = min(range_z)
            upper[z_variable] = max(range_z)
        raw_vertex_data=DSCaseVerticesFor3DSliceAndConnectivity(self._swigwrapper, 
                                                                lower._swigwrapper,
                                                                upper._swigwrapper,
                                                                x_variable,
                                                                y_variable,
                                                                z_variable)
        log_vertices = DSMatrixArrayMatrix(raw_vertex_data, 0)
        connectivity = DSMatrixArrayMatrix(raw_vertex_data, 1)
        if log_out is True:
            vertices=log_vertices
        else:
            vertices = list()
            for vertex in log_vertices:
                vertices.append([10**coordinate for coordinate in vertex])
        return [vertices, connectivity]
    
    def faces_3D_slice(self, p_vals, x_variable, y_variable, z_variable, 
                          range_x=None, range_y=None, range_z=None,
                          log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        faces = list()
        if range_x is None:
            lower[x_variable] = 1E-20
            upper[x_variable] = 1E20
        else:
            lower[x_variable] = min(range_x)
            upper[x_variable] = max(range_x)
        if range_y is None:
            lower[y_variable] = 1E-20
            upper[y_variable] = 1E20
        else:
            lower[y_variable] = min(range_y)
            upper[y_variable] = max(range_y)
        if range_y is None:
            lower[z_variable] = 1E-20
            upper[z_variable] = 1E20
        else:
            lower[z_variable] = min(range_z)
            upper[z_variable] = max(range_z)
        faces_data=DSCaseFacesFor3DSliceAndConnectivity(self._swigwrapper, 
                                                        lower._swigwrapper,
                                                        upper._swigwrapper,
                                                        x_variable,
                                                        y_variable,
                                                        z_variable)
        for i in xrange(DSMatrixArrayNumberOfMatrices(faces_data)):
            log_vertices = DSMatrixArrayMatrix(faces_data, i)
            if log_out is False:
                vertices = list()
                for vertex in log_vertices:
                    vertices.append([10**coordinate for coordinate in vertex])
            else:
                vertices = log_vertices
            faces.append(vertices)
        return faces

    def vertices_ND_slice(self, p_bounds=None, log_out=False):
        keys = []
        free_variables = []
        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        if p_bounds is None:
            p_bounds = {}
        for key in lower:           
            if key in p_bounds:
                try:
                    lower[key] = min(p_bounds[key])
                except:
                    lower[key] = p_bounds[key]
                try:
                    upper[key] = max(p_bounds[key])
                except:
                    upper[key] = p_bounds[key]
            else:
                lower[key]=1e-20
                upper[key]=1e20
            if lower[key] != upper[key]:
                keys.append(key)
                free_variables.append(lower.keys().index(key))
        vertices_data=DSCaseVerticesForNDSlice(self._swigwrapper, 
                                               lower._swigwrapper,
                                               upper._swigwrapper)
        vertices_raw = DSMatrixArrayMatrix(vertices_data, 0)
        vertices=[]
        for v in vertices_raw:
            vertices.append([v[i] for i in xrange(len(v)) if i in free_variables])
        connectivity = DSMatrixArrayMatrix(vertices_data, 1)
        return keys, vertices, connectivity
    
    def positive_roots(self, parameter_values):
        ssys = self.ssystem.remove_algebraic_constraints()
        roots = ssys.positive_roots(parameter_values)
        return roots

class CaseIntersection(object):
    
    def __init__(self, cases, name=None):
        
        setattr(self, '_name', 'Unnamed')
        setattr(self, '_cases', list())
        if name is not None:
            self._name = str(name)
        if isinstance(cases, list) is False:
            cases= [cases]
        for case in cases:
            if isinstance(case, Case) is False:
                raise TypeError, 'must be an instance of the Case class'
            self._cases.append(case)
        return
    
    def __str__(self):
        jstr = ', '
        return jstr.join([str(i) for i in self._cases])
        
    def __repr__(self):
        return 'CaseIntersection: Cases ' + str(self)
        
    def _is_valid_slice(self, p_bounds):

        lower = VariablePool(names=self._cases[0].independent_variables)
        upper = VariablePool(names=self._cases[0].independent_variables)
        for i in lower:
            lower[i] = 1E-20
            upper[i] = 1E20
        for (key,value) in p_bounds.iteritems():
            try:
                min_value,max_value = value
            except TypeError:
                min_value = value
                max_value = value
            if min_value > max_value:
                raise ValueError, 'parameter slice bounds are inverted: min is larger than max'
            lower[key] = min_value
            upper[key] = max_value
        return DSCaseIntersectionIsValidAtSlice(len(self._cases),
                                                [i._swigwrapper for i in self._cases],
                                                lower._swigwrapper,
                                                upper._swigwrapper
                                                )
    
    def is_valid(self, p_bounds=None):
        if p_bounds is not None:
            return self._is_valid_slice(p_bounds)
        return DSCaseIntersectionIsValid(len(self._cases), [i._swigwrapper for i in self._cases])
    
    def valid_parameter_set_excluding_slice(self, names):
        if isinstance(names, list) is False:
            names = [names]
        cases = self._cases
        vp=DSCaseIntersectionExceptSliceValidParameterSet(len(cases), 
                                                          [i._swigwrapper for i in cases],
                                                          len(names), names
                                                          )
        if vp is None:
            return None
        pvals = VariablePool()
        pvals.set_swigwrapper(vp)
        return pvals
        
    def vertices_1D_slice(self, p_vals, slice_variable, range_slice=None,
                          log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        if range_slice is None:
            lower[slice_variable] = 1E-20
            upper[slice_variable] = 1E20
        else:
            lower[slice_variable] = min(range_slice)
            upper[slice_variable] = max(range_slice)
        log_vertices=DSCaseIntersectionVerticesForSlice(len(self._cases),
                                                        [i._swigwrapper for i in self._cases],
                                                        lower._swigwrapper,
                                                        upper._swigwrapper,
                                                        1,
                                                        [slice_variable]
                                                        )
        vertices = list()
        for vertex in log_vertices:
            vertices.append([10**coordinate for coordinate in vertex])
        if log_out is True:
            return log_vertices
        return vertices


    def vertices_2D_slice(self, p_vals, x_variable, y_variable, range_x=None, range_y=None,
                          log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        if range_x is None:
            lower[x_variable] = 1E-20
            upper[x_variable] = 1E20
        else:
            lower[x_variable] = min(range_x)
            upper[x_variable] = max(range_x)
        if range_y is None:
            lower[y_variable] = 1E-20
            upper[y_variable] = 1E20
        else:
            lower[y_variable] = min(range_y)
            upper[y_variable] = max(range_y)
        log_vertices=DSCaseIntersectionVerticesForSlice(len(self._cases),
                                                        [i._swigwrapper for i in self._cases],
                                                        lower._swigwrapper,
                                                        upper._swigwrapper,
                                                        2,
                                                        [x_variable, y_variable]
                                                        )
        vertices = list()
        for vertex in log_vertices:
            vertices.append([10**coordinate for coordinate in vertex])
        if log_out is True:
            return log_vertices
        return vertices
    
    def faces_3D_slice(self, p_vals, x_variable, y_variable, z_variable, 
                          range_x=None, range_y=None, range_z=None,
                          log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        faces = list()
        if range_x is None:
            lower[x_variable] = 1E-20
            upper[x_variable] = 1E20
        else:
            lower[x_variable] = min(range_x)
            upper[x_variable] = max(range_x)
        if range_y is None:
            lower[y_variable] = 1E-20
            upper[y_variable] = 1E20
        else:
            lower[y_variable] = min(range_y)
            upper[y_variable] = max(range_y)
        if range_y is None:
            lower[z_variable] = 1E-20
            upper[z_variable] = 1E20
        else:
            lower[z_variable] = min(range_z)
            upper[z_variable] = max(range_z)
        faces_data=DSCaseIntersectionFacesFor3DSliceAndConnectivity(
         len(self._cases),
         [i._swigwrapper for i in self._cases], 
         lower._swigwrapper,
         upper._swigwrapper,
         x_variable,
         y_variable,
         z_variable)
        for i in xrange(DSMatrixArrayNumberOfMatrices(faces_data)):
            log_vertices = DSMatrixArrayMatrix(faces_data, i)
            if log_out is False:
                vertices = list()
                for vertex in log_vertices:
                    vertices.append([10**coordinate for coordinate in vertex])
            else:
                vertices = log_vertices
            faces.append(vertices)
        return faces

 
