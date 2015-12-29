''' Definition of the abstract model class.


'''

from __future__ import division
from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.models.ssystem import SSystem
from dspace.expressions import Expression

from math import *


class Case(Model):
    
    def __init__(self, model, swigwrapper, name=None, constraints=None, latex_symbols=None, **kwargs):
        ''' Init method for the model base class.
        
        The model object is initialized with data to construct
        a dictionary of key value pairs, where the keys represent
        names of the model's dependent variables, and the values
        are the equations in the system.
        '''
        if name == None:
            name = 'Unnamed'
        super(Case, self).__init__(model.equations,
                                   name=name,
                                   latex_symbols=latex_symbols)
        setattr(self, '_ssystem', None)
        setattr(self, '_independent_variables', None)
        self.set_swigwrapper(swigwrapper)
        if constraints is not None:
            if isinstance(constraints, list) is False:
                constraints = [constraints]
            if len(constraints) > 0:
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
            
    def __getstate__(self):
        odict = self.__dict__.copy()
        odict['_swigwrapper'] = DSSWIGDSCaseEncodedBytes(self._swigwrapper)
        del odict['_ssystem']
        del odict['_independent_variables']
        return odict
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        encoded = state['_swigwrapper']
        self.set_swigwrapper(DSSWIGDSCaseDecodeFromByteArray(encoded)) 
               
    def set_swigwrapper(self, case_swigwrapper):
        self._swigwrapper = case_swigwrapper
        Xd = VariablePool()
        Xd.set_swigwrapper(DSVariablePoolCopy(DSSSystemXd(DSCaseSSystem(case_swigwrapper))))
        for i in VariablePool():
            if i not in self.dependent_variables:
                raise NameError, 'Dependent Variables are inconsistent'
        self._dependent_variables = Xd
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
                                swigwrapper=DSCaseSSystem(case_swigwrapper),
                                latex_symbols=self._latex)
    
    @property
    def dependent_variables(self):
        return self._dependent_variables.keys()
        
    @property
    def equations(self):
        eqs = DSCaseEquations(self._swigwrapper)
        equations = list()
        for i in xrange(0, DSCaseNumberOfEquations(self._swigwrapper)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equations.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        return Equations(equations, latex_symbols=self._latex)
    
    @property
    def ssystem(self):
        return self._ssystem
    
    @property
    def independent_variables(self):
        return self._independent_variables.keys()
    
    def _valid_parameter_set_bounded(self, p_bounds, optimize=None, minimize=True):
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
        if optimize is not None:
            variablepool = DSCaseValidParameterSetAtSliceByOptimizingFunction(self._swigwrapper, 
                                                                              lower._swigwrapper,
                                                                              upper._swigwrapper,
                                                                              optimize,
                                                                              minimize)
        else:
            variablepool = DSCaseValidParameterSetAtSlice(self._swigwrapper, 
                                                          lower._swigwrapper,
                                                          upper._swigwrapper)
        pvals = VariablePool()
        pvals.set_swigwrapper(variablepool)
        return pvals
    
    def valid_parameter_set(self, p_bounds=None, optimize=None, minimize=True, strict=True):
        if p_bounds is not None:
            pvals = self._valid_parameter_set_bounded(p_bounds, 
                                                      optimize=optimize,
                                                      minimize=minimize)
            return pvals
        if optimize is not None:
            variablepool = DSCaseValidParameterSetByOptimizingFunction(self._swigwrapper, 
                                                                       optimize,
                                                                       minimize)
        else:
            variablepool = DSCaseValidParameterSet(self._swigwrapper)
        pvals = VariablePool()
        pvals.set_swigwrapper(variablepool)
        for i in pvals:
            if pvals[i] == 0.0 or pvals[i] == float('inf'):
                return self._valid_parameter_set_bounded({i:[1e-20, 1e20] for i in self.independent_variables}, 
                                                         optimize=optimize,
                                                         minimize=minimize)
        return pvals
    
    def valid_interior_parameter_set(self, p_bounds=None, distance=50, **kwargs):
        pvals = self.valid_parameter_set(p_bounds=p_bounds, **kwargs)
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

    def consistent_parameter_and_state(self):
        variablepool = DSCaseConsistentParameterAndStateSet(self._swigwrapper)
        pstate = VariablePool()
        pstate.set_swigwrapper(variablepool)
        return pstate
    
    def valid_parameter_and_state(self):
        variablepool = DSCaseValidParameterAndStateSet(self._swigwrapper)
        pstate = VariablePool()
        pstate.set_swigwrapper(variablepool)
        return pstate
          
    @property
    def case_number(self):
        return DSCaseIdentifier(self._swigwrapper)

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
        if eqs_expr is None:
            return
        for i in xrange(0, DSCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(conditions, latex_symbols=self._latex)
        
    @property
    def conditions_log(self):
        conditions = list()
        eqs_expr = DSCaseLogarithmicConditions(self._swigwrapper)
        if eqs_expr is None:
            return
        for i in xrange(0, DSCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(conditions, latex_symbols=self._latex)
    
    @property
    def boundaries(self):
        boundaries = list()
        eqs_expr = DSCaseBoundaries(self._swigwrapper)
        if eqs_expr is None:
            return
        for i in xrange(0, DSCaseNumberOfBoundaries(self._swigwrapper)):
            boundaries.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(boundaries, latex_symbols=self._latex)
    
    @property
    def boundaries_log(self):
        boundaries = list()
        eqs_expr = DSCaseLogarithmicBoundaries(self._swigwrapper)
        if eqs_expr is None:
            return
        for i in xrange(0, DSCaseNumberOfBoundaries(self._swigwrapper)):
            boundaries.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(boundaries, latex_symbols=self._latex)
    
    @property
    def is_cyclical(self):
        return False
    
    def _is_valid_slice(self, p_bounds, strict=True):

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
                                    upper._swigwrapper,
                                    strict)
    
    def steady_state(self, parameter_values):
        return self.ssystem.steady_state(parameter_values)
        
    def steady_state_flux(self, parameter_values):
        return self.ssystem.steady_state_flux(parameter_values)
        
    def steady_state_function(self, function, parameter_values):
        return self.ssystem.steady_state_function(function, parameter_values)
    
    def is_consistent(self, point=None):
        
        if point is not None:
            p_vals_d = VariablePool(names=self.dependent_variables)
            p_vals_i = VariablePool(names=self.independent_variables)
            for i in p_vals_d:
                p_vals_d[i] = point[i]
            for i in p_vals_i:
                p_vals_i[i] = point[i]
            return DSCaseIsValidInStateSpaceAtPoint(self._swigwrapper,
                                                    p_vals_d._swigwrapper,
                                                    p_vals_i._swigwrapper)
            ## return self._is_valid_slice(p_bounds)
            #do something
        return DSCaseConditionsAreValid(self._swigwrapper)
 
    def is_valid(self, p_bounds=None, strict=True):
        
        if p_bounds is not None:
            return self._is_valid_slice(p_bounds, strict=strict)
            #do something
        return DSCaseIsValid(self._swigwrapper, strict)
        
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
        dependent.update(v_bounds)
        if len(ssys.auxiliary_variables) > 0:
            aux = ssys.value_for_auxiliary_variables(Xd_t, independent)
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
            
    def vertices_2D_slice(self, p_vals, x_variable, y_variable, range_x=None,
                          range_y=None, log_out=False, vtype='numerical'):
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
        if vtype.lower() in ['numerical', 'both'] or vtype.lower() in ['n', 'b']:
            log_vertices=DSCaseVerticesFor2DSlice(self._swigwrapper, 
                                                  lower._swigwrapper,
                                                  upper._swigwrapper,
                                                  x_variable,
                                                  y_variable)
            vertices = list()
            for vertex in log_vertices:
                vertices.append([10**coordinate for coordinate in vertex])
            if log_out is True:
                vertices=log_vertices
        if vtype.lower() == 'analytical' or vtype.lower() == 'a':
            vertices=self._vertex_equations_2D_slice(p_vals, x_variable, y_variable,
                                            range_x, range_y, log_out)
        if vtype.lower() == 'both' or vtype.lower() == 'b':
            
            eqs=self._vertex_equations_2D_slice(p_vals, x_variable, y_variable,
                                                range_x, range_y, log_out)
            vertices=[(vertices[i],eqs[i]) for i in xrange(len(vertices))]
        return vertices

    def _vertex_equations_2D_slice(self, p_vals, x_variable, y_variable, range_x, range_y,
                                  log_out):
        lower = p_vals.copy()
        upper = p_vals.copy()
        lower[x_variable] = min(range_x)
        upper[x_variable] = max(range_x)
        lower[y_variable] = min(range_y)
        upper[y_variable] = max(range_y)
        stack = DSCaseVertexEquationsFor2DSlice(self._swigwrapper, 
                                                lower._swigwrapper,
                                                upper._swigwrapper,
                                                x_variable,
                                                y_variable,
                                                log_out)
        vertices = []
        for i in xrange(DSStackCount(stack)):
            expressions = DSExpressionArrayFromVoid(DSStackPop(stack))
            expr0 = Expression(None)
            expr0._swigwrapper = DSExpressionAtIndexOfExpressionArray(expressions, 0)
            expr1 = Expression(None)
            expr1._swigwrapper = DSExpressionAtIndexOfExpressionArray(expressions, 1)
            vertices.append(Equations([str(expr0), str(expr1)], latex_symbols=self._latex))
            DSSecureFree(expressions)
        DSStackFree(stack)
        vertices.reverse()
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
    
    def eigen_spaces(self):
        ds_swig = DSCaseEigenSubspaces(self._swigwrapper)
        eqs = DSDesignSpaceEquations(ds_swig)
        equations = list()
        for i in xrange(0, DSDesignSpaceNumberOfEquations(ds_swig)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equations.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        import dspace.models.designspace as designspace
        ds = designspace.DesignSpace(Equations(equations),
                                     latex_symbols=self._latex,
                                     swigwrapper=ds_swig)
        return ds


class CaseIntersection(Case):
    
    def __init__(self, cases, name=None, constraints=None, latex_symbols=None):
        
        
        setattr(self, '_name', 'Unnamed')
        setattr(self, '_cases', list())
        setattr(self, '_latex', {})
        if isinstance(cases, list) is False:
            cases= [cases]
        for case in cases:
            if isinstance(case, Case) is False:
                raise TypeError, 'must be an instance of the Case class'
        self._latex = case._latex
        new_cases = []
        if constraints is not None:
            if isinstance(constraints, list) is False:
                constraints = [constraints]
        if constraints is not None:
            cases_swig = [DSCaseCopy(i._swigwrapper) for i in cases]
            for i in range(len(cases_swig)):
                DSCaseAddConstraints(cases_swig[i], constraints, len(constraints))
                new_cases.append(Case(cases[i], cases_swig[i], name=cases[i].name))
        else:
            new_cases = cases
        swigwrapper = DSSWIGPseudoCaseFromIntersectionOfCases(len(cases), [i._swigwrapper for i in new_cases])
        super(CaseIntersection, self).__init__(cases[0],
                                               swigwrapper,
                                               name=name,
                                               latex_symbols=latex_symbols)
        self._cases = new_cases
        return

    def __del__(self):
        ''' 
        '''
        if self._swigwrapper is not None:
            DSCaseFree(self._swigwrapper)
            
    def set_swigwrapper(self, case_swigwrapper):
        self._swigwrapper = case_swigwrapper
        Xd = VariablePool()
        Xd.set_swigwrapper(DSVariablePoolCopy(DSCaseXd(case_swigwrapper)))
        for i in VariablePool():
            if i not in self.dependent_variables:
                raise NameError, 'Dependent Variables are inconsistent'
        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSCaseXi(case_swigwrapper)))
        self._independent_variables = Xi
    
    def __str__(self):
        jstr = ', '
        return jstr.join([str(i) for i in self._cases])
    ##     
    def __repr__(self):
        return 'CaseIntersection: Cases ' + str(self)
        
class CaseColocalization(CaseIntersection):
    
    def __init__(self, cases, slice_variables, name=None, constraints=None, latex_symbols=None):
        
        
        setattr(self, '_name', 'Unnamed')
        setattr(self, '_cases', list())
        setattr(self, '_latex', {})
        slice_constraints = []
        case_constraints = []
        if isinstance(cases, list) is False:
            cases= [cases]
        if isinstance(slice_variables, list) is False:
            slice_variables = [slice_variables]
        for case in cases:
            if isinstance(case, Case) is False:
                raise TypeError, 'must be an instance of the Case class'
        self._latex = case._latex
        new_cases = []
        if constraints is not None:
            if isinstance(constraints, list) is False:
                constraints = [constraints]
            for i in constraints:
                if '$' in i:
                    slice_constraints.append(i)
                else:
                    case_constraints.append(i)
        if len(case_constraints) > 0:
            cases_swig = [DSCaseCopy(i._swigwrapper) for i in cases]
            for i in range(len(cases_swig)):
                DSCaseAddConstraints(cases_swig[i], case_constraints, len(case_constraints))
                new_cases.append(Case(cases[i], cases_swig[i], name=cases[i].name))
        else:
            new_cases = cases
        swigwrapper = DSSWIGPseudoCaseFromIntersectionOfCasesExcludingSlice(len(cases),
                                                                            [i._swigwrapper for i in new_cases],
                                                                            len(slice_variables),
                                                                            slice_variables)
        if len(slice_constraints) > 0:
            DSCaseAddConstraints(swigwrapper, slice_constraints, len(slice_constraints))
        super(CaseIntersection, self).__init__(cases[0],
                                                 swigwrapper,
                                                 name=name,
                                                 latex_symbols=latex_symbols)
        self._cases = new_cases
        self._slice_variables = slice_variables
        return

    def __del__(self):
        ''' 
        '''
        if self._swigwrapper is not None:
            DSCaseFree(self._swigwrapper)
            
    def set_swigwrapper(self, case_swigwrapper):
        self._swigwrapper = case_swigwrapper
        Xd = VariablePool()
        Xd.set_swigwrapper(DSVariablePoolCopy(DSCaseXd(case_swigwrapper)))
        for i in VariablePool():
            if i not in self.dependent_variables:
                raise NameError, 'Dependent Variables are inconsistent'
        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSCaseXi(case_swigwrapper)))
        self._independent_variables = Xi
    
    def __str__(self):
        jstr = ', '
        return jstr.join([str(i) for i in self._cases])
        
    def __repr__(self):
        return 'CaseColocalization: Cases ' + str(self)
        
    def valid_parameter_set(self, p_bounds=None, optimize=None, minimize=True, project=True):
        psetHD = super(CaseIntersection, self).valid_parameter_set(
                                          p_bounds=p_bounds,
                                          optimize=optimize,
                                          minimize=minimize)
        if len(psetHD) == 0:
            return psetHD            
        if project is False:
            return psetHD
        p_sets = dict()
        index = 0
        for i in self._cases:
            pvals = VariablePool(names=i.independent_variables)
            for j in pvals:
                if j in self._slice_variables:
                    pvals[j] = psetHD['$'+j+'_'+str(index)]
                else:
                    pvals[j] = psetHD[j]
            p_sets[str(i)] = pvals
            index += 1
        return p_sets
        
    def valid_interior_parameter_set(self, distance=50, p_bounds=None, project=True, **kwargs):
        psetHD = super(CaseIntersection, self).valid_interior_parameter_set(
                                          p_bounds=p_bounds,
                                          distance=distance,
                                          project=False, 
                                          **kwargs)
        if len(psetHD) == 0:
            return psetHD
        if project is False:
            return psetHD
        p_sets = dict()
        index = 0
        for i in self._cases:
            pvals = VariablePool(names=i.independent_variables)
            for j in pvals:
                if j in self._slice_variables:
                    pvals[j] = psetHD['$'+j+'_'+str(index)]
                else:
                    pvals[j] = psetHD[j]
            p_sets[str(i)] = pvals
            index += 1
        return p_sets
