''' Definition of the abstract model class.


'''

from __future__ import division
from itertools import chain
from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.models.ssystem import SSystem
from dspace.expressions import Expression
import numpy as np
import itertools
from math import *
import random


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
                raise NameError('Dependent Variables are inconsistent')
        self._dependent_variables = Xd
        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSSSystemXi(DSCaseSSystem(case_swigwrapper))))
        self._independent_variables = Xi
        eqs = list()
        eqs_expr = DSSSystemEquations(DSCaseSSystem(case_swigwrapper))

        #for i in xrange(0, len(super(Case, self).equations)): #revert this! Eventually use DSVariablePoolNumberOfVariables(DSCaseXd(aCase))
        for i in range(0, DSVariablePoolNumberOfVariables(DSCaseXd(case_swigwrapper))):
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
        for i in range(0, DSCaseNumberOfEquations(self._swigwrapper)):
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
    
    def valid_parameter_set(self, p_bounds=None, optimize=None, minimize=True, strict=True, **kwargs):
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
            if 'shared_boundaries' in kwargs:
                variablepool = DSCaseSharedBoundariesValidParameterSet(self._swigwrapper)
            else:
                variablepool = DSCaseValidParameterSet(self._swigwrapper)
        pvals = VariablePool()
        pvals.set_swigwrapper(variablepool)
        for i in pvals:
            # if pvals[i] == 0.0 or pvals[i] == float('inf'): #jason original
            if pvals[i] == 0.0 or pvals[i] > 1e20 or pvals[i] < 1e-20:  # Miguel edit
                return self._valid_parameter_set_bounded({i: [1e-20, 1e20] for i in self.independent_variables},
                                                         optimize=optimize,
                                                         minimize=minimize)

        return pvals
    
    def valid_interior_parameter_set(self, p_bounds=None, distance=50, **kwargs):

        pvals = self.valid_parameter_set(p_bounds=p_bounds, **kwargs)

        for j in range(0, 2):
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
                    if len(slice_range) > 1:
                        pvals[i] = 10**(slice_range[0] + (slice_range[1] - slice_range[0])/2)
                    else:
                        pvals[i] = slice_range[0]

        return pvals

    def valid_interior_parameter_bounding_box(self, p_bounds=None, log_out=False):
        tolerances = self.bounding_box(p_bounds=p_bounds, log_out=True)
        pvals = VariablePool(names=self.independent_variables)
        for key, values in tolerances.items():
            pvals[key] = (values[0] + values[1])/2
            if log_out is False:
                pvals[key] = 10**pvals[key]
        return pvals

    def valid_interior_parameter_set_vertex_enumeration(self, p_bounds=None, log_out=False):
        pvals = VariablePool(names=self.independent_variables)
        lb = pvals.copy()
        ub = pvals.copy()
        for key in lb.keys():
            if p_bounds is None:
                lb[key] = 1e-20
                ub[key] = 1e20
            else:
                lb[key] = p_bounds[key][0]
                ub[key] = p_bounds[key][1]
        maxVertices = 0
        limitVertices = False
        op_swig = DSCaseCentroid_qhull(self._swigwrapper,
                                       lb._swigwrapper,
                                       ub._swigwrapper,
                                       maxVertices,
                                       limitVertices)

        operating_point = VariablePool(names=self.independent_variables)
        for var in operating_point.keys():
            operating_point[var] = DSVariablePoolValueForVariableWithName(op_swig, var)
        DSVariablePoolFree(op_swig)

        for key, values in operating_point.items():
            if log_out is True:
                pvals[key] = operating_point[key]
            else:
                pvals[key] = 10**operating_point[key]

        return pvals

    def valid_interior_parameter_set_vertex_enumeration_average_centroid(self, p_bounds=None, log_out=False):
        pvals = VariablePool(names=self.independent_variables)
        lb = pvals.copy()
        ub = pvals.copy()
        for key in lb.keys():
            if p_bounds is None:
                lb[key] = 1e-20
                ub[key] = 1e20
            else:
                lb[key] = p_bounds[key][0]
                ub[key] = p_bounds[key][1]
        maxVertices = 0
        limitVertices = False
        volume, nr_vertices, vertices_matrix, operating_point = self.volume_lrs(lb,
                                                                                 ub,
                                                                                 maxVertices,
                                                                                 limitVertices,
                                                                                 return_vertices_matrix=True)
        for key, values in operating_point.items():
            if log_out is True:
                pvals[key] = operating_point[key]
            else:
                pvals[key] = 10**operating_point[key]

        return pvals

    def valid_interior_parameter_geometric_mean_t_bb(self, p_bounds=None, log_out=False, distance=50, **kwargs):
        pvals_tolerance = self.valid_interior_parameter_set(p_bounds=p_bounds,
                                                            log_out=log_out,
                                                            distance=distance,
                                                            **kwargs)
        pvals_bb = self.valid_interior_parameter_bounding_box(p_bounds=p_bounds,
                                                              log_out=log_out)

        pvals_gmean = pvals_tolerance.copy()
        for key in pvals_gmean:
            pvals_gmean[key] = (pvals_tolerance[key]*pvals_bb[key])**0.5
        return pvals_gmean

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

    def volume(self, ignore_unbounded=True, log_coordinate=False, shared_boundaries=False, method='Tolerances',
               **kwargs):

        neg_vol = False

        # 06 April 2021. Including p_bounds in the identification of pvals guarantees that the operating point of
        # the system for which tolerances are calculated is located within the boundaries being analyzed.

        if 'lowerBounds' in kwargs.keys() and 'upperBounds' in kwargs.keys():
            lowerBounds = kwargs['lowerBounds']
            upperBounds = kwargs['upperBounds']
            p_bounds = dict((i, [lowerBounds, upperBounds]) for i in self.independent_variables)
            asymmetrical = False
        elif 'p_bounds' in kwargs.keys():
            p_bounds = kwargs['p_bounds']
            asymmetrical = True
        else:
            raise NameError('Individual bounds (upperBounds and lowerBounds) or p_bounds should be provided')

        pvals = self.valid_interior_parameter_set(p_bounds=p_bounds,
                                                  shared_boundaries=shared_boundaries)

        if method == 'Tolerances':
            try:
                tolerances = self.measure_tolerance(pvals,
                                                    log_out=log_coordinate,
                                                    shared_boundaries=shared_boundaries)
            except:
                try:
                    p_bounds = dict((i, [10**-20, 10**20]) for i in pvals)
                    pvals = self.valid_interior_parameter_set(p_bounds=p_bounds, shared_boundaries=shared_boundaries)
                    tolerances = self.measure_tolerance(pvals, log_out=log_coordinate,
                                                        shared_boundaries=shared_boundaries)
                except:
                    print("Calculation of tolerances for case {} failed. Assuming a volume of 1".format(self.case_number))
                    volume = 1
                    return volume

        elif method == 'Bounding Box':
            tolerances = self.bounding_box(log_out=log_coordinate)

        elif method == 'Geometric Mean T. & BB.':

            volume_tol = self.volume(ignore_unbounded=ignore_unbounded,
                                 log_coordinate=log_coordinate,
                                 shared_boundaries=shared_boundaries,
                                 method='Tolerances',
                                 **kwargs)
            volume_bb = self.volume(ignore_unbounded=ignore_unbounded,
                                 log_coordinate=log_coordinate,
                                 shared_boundaries=shared_boundaries,
                                 method='Bounding Box',
                                 **kwargs)

            vol = (volume_tol * volume_bb)**0.5
            return vol
        else:                                           #lrs library
            maxVertices = kwargs['maxVertices']
            limitVertices = kwargs['limitVertices']
            lb = pvals.copy()
            ub = pvals.copy()
            for key in lb.keys():
                if asymmetrical is False:
                    lb[key] = lowerBounds
                    ub[key] = upperBounds
                else:
                    lb[key] = min(p_bounds[key])
                    ub[key] = max(p_bounds[key])
            volume, vertices = self.volume_lrs(lb, ub, maxVertices, limitVertices)
            return volume

        volume = 1
        for xi in sorted(pvals.keys()):

            if asymmetrical is True:
                lowerBounds = min(p_bounds[xi])
                upperBounds = max(p_bounds[xi])

            lower_th = lowerBounds if log_coordinate is False else log10(lowerBounds)  # 1e-15
            upper_th = upperBounds if log_coordinate is False else log10(upperBounds)   # 1e15
            lower, upper = tolerances[xi]

            if method == 'Tolerances':
                lower = lower*pvals[xi] if log_coordinate is False else lower + log10(pvals[xi])
                upper = upper*pvals[xi] if log_coordinate is False else upper + log10(pvals[xi])

            lower_value = lower if lower > lower_th else lower_th
            upper_value = upper if upper < upper_th else upper_th

            ratio = upper_value/lower_value if log_coordinate is False else upper_value - lower_value

            if ratio < 0.0:
                neg_vol = True
            if shared_boundaries is True and ratio == 0.0:
                ratio = 1
                print("Case {}: Changing ratio for the variable {} from zero to {}".format(self.case_number, xi, ratio))
            if ignore_unbounded is True:
                if lower > lower_th and upper < upper_th:
                    volume = volume*ratio
            else:
                volume = volume * ratio
        volume = round_sig(volume) if neg_vol is False and volume != 0.0 else 0.0
        return volume

    def volume_geometric_mean(self, ignore_unbounded=True, log_coordinate=False):
        pvals = self.valid_interior_parameter_set()
        tolerances = self.measure_tolerance(pvals, log_out=log_coordinate)
        parameter_count = 0
        volume = 1
        for xi in sorted(pvals.keys()):
            lower_th = 1e-15 if log_coordinate is False else -15
            upper_th = 1e15 if log_coordinate is False else 15
            lower, upper = tolerances[xi]
            if log_coordinate is False:
                lower_value = lower if lower > lower_th else 1e-3
                upper_value = upper if upper < upper_th else 1e3
            else:
                lower_value = lower if lower > lower_th else -3
                upper_value = upper if upper < upper_th else 3
            ratio = upper_value / lower_value if log_coordinate is False else upper_value - lower_value
            if ignore_unbounded is True:
                if lower > lower_th and upper < upper_th:
                    volume = volume*ratio
                    parameter_count += 1
            else:
                volume = volume * ratio
                parameter_count += 1
        if parameter_count == 0:
            parameter_count = 1
        geometric_mean = round_sig(volume**(1.0/parameter_count))
        return geometric_mean

    def volume_lrs(self, lowerBounds, upperBounds, maxVertices, limitVertices, return_vertices_matrix=False):
        vol_structure = DSCaseVolume_lrs(self._swigwrapper, lowerBounds._swigwrapper,
                                upperBounds._swigwrapper, maxVertices, limitVertices, return_vertices_matrix)
        volume = DSCaseVolumeGetVolume(vol_structure)
        nr_vertices = DSCaseVolumeGetVertices(vol_structure)

        if return_vertices_matrix is False:
            return volume, nr_vertices
        else:
            vertices_matrix = DSCaseVolumeGetVerticesMatrix(vol_structure)
            op_swig = DSCaseVolumeGetOperatingPoint2D(vol_structure)
            operating_point = lowerBounds.copy()
            for var in operating_point:
                operating_point.update({var: DSVariablePoolValueForVariableWithName(op_swig, var)})
            return volume, nr_vertices, vertices_matrix, operating_point

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
            return None
        for i in range(0, DSCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(conditions, latex_symbols=self._latex)
        
    @property
    def conditions_log(self):
        conditions = list()
        eqs_expr = DSCaseLogarithmicConditions(self._swigwrapper)
        if eqs_expr is None:
            return None
        for i in range(0, DSCaseNumberOfConditions(self._swigwrapper)):
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
        for i in range(0, DSCaseNumberOfBoundaries(self._swigwrapper)):
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
        for i in range(0, DSCaseNumberOfBoundaries(self._swigwrapper)):
            boundaries.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(boundaries, latex_symbols=self._latex)
    
    @property
    def is_cyclical(self):
        return False

    @property
    def is_unstable(self):
        return DSSSystemIsUnstable(DSCaseSSystem(self._swigwrapper))

    @property
    def is_false_blowing(self):
        return DSSSystemIsFalseBlowing(DSCaseSSystem(self._swigwrapper))

    def adjust_codominant_stoichiometry(self):
        DSCaseAdjustStoichiometryOfCodominantCase(self._swigwrapper)

    @property
    def should_adjust_codominant_stoichiometry(self):
        return DSSSystemAdjustCodominantStoichiometry(DSCaseSSystem(self._swigwrapper))

    @property
    def number_of_mass_balances(self):
        return DSCaseNumberOfMassBalances(self._swigwrapper)

    def _is_valid_slice(self, p_bounds, strict=True):

        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        for i in lower:
            lower[i] = 1E-20
            upper[i] = 1E20
        for (key, value) in p_bounds.items():
            try:
                min_value, max_value = value
            except TypeError:
                min_value = value
                max_value = value
            if min_value > max_value:
                raise ValueError('parameter slice bounds are inverted: min is larger than max')
            lower[key] = min_value
            upper[key] = max_value
        return DSCaseIsValidAtSlice(self._swigwrapper,
                                    lower._swigwrapper,
                                    upper._swigwrapper,
                                    strict)

    def dominant_fin_at_index(self, i):
        return DSCaseDominantFinAtIndex(self._swigwrapper, i)

    def dominant_fout_at_index(self, i):
        return DSCaseDominantFoutAtIndex(self._swigwrapper, i)
    
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
            for key, value in p_bounds.items():
                if key not in lower:
                    continue;
                try:
                    min_value,max_value = value
                except TypeError:
                    min_value = value
                    max_value = value
                    boundkeys.remove(key)
                if min_value > max_value:
                    raise ValueError('Min cannot be larger than max')
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
            
    def measure_tolerance(self, pvals, log_out=False, shared_boundaries=False):
        tolerances = {}
        for key in self.independent_variables:
            tolerances[key] = self.vertices_1D_slice(pvals, key, log_out=log_out)
            dimension = len(tolerances[key])
            if shared_boundaries is True and dimension == 1:
                tolerances[key] = [tolerances[key][0], tolerances[key][0]]
            if log_out is False:
                tolerances[key] = (tolerances[key][0][0]/pvals[key], tolerances[key][1][0]/pvals[key])
            if log_out is True:
                tolerances[key] = (tolerances[key][0][0]-log10(pvals[key]), tolerances[key][1][0]-log10(pvals[key]))
        return tolerances

    def double_value_boundaries_at_point(self, pvals, log_out=False):
        boundaries = list(chain(*DSCaseDoubleValueBoundariesAtPointSortXi(self._swigwrapper, pvals._swigwrapper)))
        if log_out is True:
            return boundaries
        else:
            boundaries_cartesian = [pow(10, i) for i in boundaries]
            return boundaries_cartesian
              
    def vertices_1D_slice(self, p_vals, slice_variable, range_slice=None, log_out=False):
        lower = p_vals.copy()
        upper = p_vals.copy()
        if range_slice is None:
            lower[slice_variable] = 1E-20 # original values were 1e-20 // modified values were 1e-25
            upper[slice_variable] = 1E20  # original values were 1e20 // modified values were 1e20
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
            V = list(zip(*V))[0]
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
                vertices = log_vertices
        if vtype.lower() == 'analytical' or vtype.lower() == 'a':
            vertices = self._vertex_equations_2D_slice(p_vals, x_variable, y_variable,
                                            range_x, range_y, log_out)
        if vtype.lower() == 'both' or vtype.lower() == 'b':
            
            eqs = self._vertex_equations_2D_slice(p_vals, x_variable, y_variable,
                                                range_x, range_y, log_out)
            vertices = [(vertices[i],eqs[i]) for i in range(len(vertices))]
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
        for i in range(DSStackCount(stack)):
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
        for i in range(DSMatrixArrayNumberOfMatrices(faces_data)):
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
            vertices.append([v[i] for i in range(len(v)) if i in free_variables])
        connectivity = DSMatrixArrayMatrix(vertices_data, 1)
        return keys, vertices, connectivity
    
    def positive_roots(self, parameter_values):
        ssys = self.ssystem.remove_algebraic_constraints()
        roots = ssys.positive_roots(parameter_values)
        return roots

    def positive_roots_numpy(self, parameter_values):
        ssys = self.ssystem.remove_algebraic_constraints()
        # try:
        roots = ssys.positive_roots_numpy(parameter_values)
        # except:
        #     print("Case {} caused the routine positive_roots_numpy to crash".format(self.case_number))
        #     return 0
        return roots

    def has_complex_conjugates(self, parameter_values):
        ssys = self.ssystem.remove_algebraic_constraints()
        try:
            complex_conjugates = ssys.has_complex_conjugates(parameter_values)
        except:
            print("Case {} caused the routine has_complex_conjugates to crash".format(self.case_number))
            return False
        return complex_conjugates
    
    def eigen_spaces(self):
        ds_swig = DSCaseEigenSubspaces(self._swigwrapper)
        eqs = DSDesignSpaceEquations(ds_swig)
        equations = list()
        for i in range(0, DSDesignSpaceNumberOfEquations(ds_swig)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equations.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        import dspace.models.designspace as designspace
        ds = designspace.DesignSpace(Equations(equations),
                                     latex_symbols=self._latex,
                                     swigwrapper=ds_swig)
        return ds

    def shared_boundaries_indices(self, case2, intersecting=False):
        return DSCaseSharedBoundaries(self._swigwrapper, case2._swigwrapper, intersecting)

    def share_boundaries_with(self, case2, intersecting=False):
        return DSCaseHasSharedBoundaries(self._swigwrapper, case2._swigwrapper, intersecting)

    def shared_boundaries_number_of_vertices(self, case2, lowerBounds, upperBounds, maxVertices, limitVertices):
        return DSCaseSharedBoundariesNumberOfVertices(self._swigwrapper, case2._swigwrapper,
                                                      lowerBounds._swigwrapper,
                                                      upperBounds._swigwrapper,
                                                      maxVertices,
                                                      limitVertices)

    def shared_boundaries_is_valid(self, case2):
        return DSCasesSharedBoundariesIsValid(self._swigwrapper, case2._swigwrapper)

    def calculate_distance_to_case(self, case2):
        p1 = self.valid_interior_parameter_set()
        p2 = case2.valid_interior_parameter_set()
        return DSVariablePoolDistanceToPool(p1._swigwrapper, p2._swigwrapper)

    def dimension(self, lowerBounds, upperBounds):
        return DSCaseDimension(self._swigwrapper, lowerBounds._swigwrapper, upperBounds._swigwrapper)

    def neighbors_signature(self, controller):
        neighbors_vector = DSCaseGetSignatureNeighbors(self._swigwrapper, controller._swigwrapper)
        neighbors_list = [str(DSUIntegerVectorValueAtIndex(neighbors_vector, i)) for i in range(DSUIntegerVectorDimension(neighbors_vector))]
        return neighbors_list

    def sample_valid_points(self, p_bounds=None, nr_points=10):
        # 1. Generate a list containing bounding box in each dimension.
        box = self.bounding_box(p_bounds=p_bounds, log_out=True)
        variables = self.independent_variables

        # 2. generate a list of coordinates for each dimension to sample the space
        coordinates_seed = []
        for key in variables:
            coordinates_seed.append(np.linspace(box[key][0], box[key][1], nr_points))

        # 3. Use the itertools.product function to generate list of points
        points = itertools.product(*coordinates_seed)

        valid_points = []
        # 4. Check for validity
        for point in points:
            validity, p = self.generate_pool_validity(point)
            if validity is True:
                valid_points.append(p)

        return valid_points

    def generate_pool_validity(self, point):
        variables = self.independent_variables
        pvals = VariablePool(names=variables)
        for index, variable in enumerate(variables):
            pvals[variable] = 10 ** point[index]
        validity = self.is_valid(p_bounds=pvals, strict=True)
        return validity, pvals

    def mutation_rate_to_phenotype_grid(self, case2, identity,
                                        lamb=1.0,
                                        delta=2.0,
                                        p_bounds=None,
                                        nr_points=0,
                                        average_method='Average'):

        # 1Get the bounding box and sample nr_points. Let us create a function for that. Do this for both case1
        # and case2.
        valid_p1 = self.sample_valid_points(p_bounds=p_bounds, nr_points=nr_points)
        valid_p2 = case2.sample_valid_points(p_bounds=p_bounds, nr_points=nr_points)

        #  For each point in valid_p1, calculate mutation rate for valid_p2.
        mutation_rate = 0
        total_points = 0

        if average_method == 'Average':
                for p1 in valid_p1:
                    # validity, p1 = self.generate_pool_validity(point1)
                    # if validity is False:
                    #     continue
                    for p2 in valid_p2:
                        # validity, p2 = case2.generate_pool_validity(point2)
                        # if validity is False:
                        #     continue
                        mutation = DSPopDynamicsMutationRateForTransition(p1._swigwrapper,
                                                                          p2._swigwrapper,
                                                                          lamb,
                                                                          delta,
                                                                          identity._swigwrapper)
                        mutation_rate += mutation
                        total_points += 1
                mutation_rate = mutation_rate/total_points
        else:
                for p1 in valid_p1:
                    # validity, p1 = self.generate_pool_validity(point1)
                    # if validity is False:
                    #     continue
                    for p2 in valid_p2:
                        # validity, p2 = case2.generate_pool_validity(point2)
                        # if validity is False:
                        #     continue
                        mutation = DSPopDynamicsMutationRateForTransition(p1._swigwrapper,
                                                                          p2._swigwrapper,
                                                                          lamb,
                                                                          delta,
                                                                          identity._swigwrapper)
                        mutation_rate += log10(mutation)
                        total_points += 1
                mutation_rate = 10**(mutation_rate/total_points)

        print("{} total mutation rates analyzed".format(total_points))
        print("The mutation rate is: ", mutation_rate)
        return mutation_rate

    def mutation_rate_to_phenotype(self, case2, identity,
                                   lamb=1.0,
                                   delta=1.5,
                                   method='Tolerances',
                                   p_bounds=None,
                                   nr_points=0,
                                   average_method='Average'):

        if identity is None:
            raise ValueError("Please set the identity of each parameter in the Edit Parameters tab")

        if method == 'Tolerances':
            p1 = self.valid_interior_parameter_set(p_bounds=p_bounds)
            p2 = case2.valid_interior_parameter_set(p_bounds=p_bounds)
        elif method == 'Bounding Box':
            p1 = self.valid_interior_parameter_bounding_box(p_bounds=p_bounds)
            p2 = case2.valid_interior_parameter_bounding_box(p_bounds=p_bounds)
        elif method == 'Vertex Enumeration':
            p1 = self.valid_interior_parameter_set_vertex_enumeration(p_bounds=p_bounds)
            p2 = case2.valid_interior_parameter_set_vertex_enumeration(p_bounds=p_bounds)
        elif method == 'Geometric Mean T. & BB.':
            p1 = self.valid_interior_parameter_geometric_mean_t_bb(p_bounds=p_bounds)
            p2 = case2.valid_interior_parameter_geometric_mean_t_bb(p_bounds=p_bounds)
        elif method == 'Grid':
            return self.mutation_rate_to_phenotype_grid(case2, identity,
                                                        lamb=lamb,
                                                        delta=delta,
                                                        p_bounds=p_bounds,
                                                        nr_points=nr_points,
                                                        average_method=average_method)

        mutation = DSPopDynamicsMutationRateForTransition(p1._swigwrapper,
                                                          p2._swigwrapper,
                                                          lamb,
                                                          delta,
                                                          identity._swigwrapper)

        return mutation

    def mutation_key_to_case(self, case2, symbol=None):
        key = 'k'

        try:
            if int(self.case_number) < 10:
                key += '0' + str(self.case_number)
            else:
                key += str(self.case_number)
        except:
            key += str(self.case_number)

        if symbol is not None:
            key += str(symbol)

        try:
            if int(case2.case_number) < 10:
                key += '0' + str(case2.case_number)
            else:
                key += str(case2.case_number)
        except:
            key += str(case2.case_number)

        return key


class CaseIntersection(Case):
    
    def __init__(self, cases, name=None, constraints=None, latex_symbols=None):
        
        
        setattr(self, '_name', 'Unnamed')
        setattr(self, '_cases', list())
        setattr(self, '_latex', {})
        if isinstance(cases, list) is False:
            cases= [cases]
        for case in cases:
            if isinstance(case, Case) is False:
                raise TypeError('must be an instance of the Case class')
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
                raise NameError('Dependent Variables are inconsistent')
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
                raise TypeError('must be an instance of the Case class')
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
                raise NameError('Dependent Variables are inconsistent')
        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSCaseXi(case_swigwrapper)))
        self._independent_variables = Xi
    
    def __str__(self):
        jstr = ', '
        return jstr.join([str(i) for i in self._cases])
        
    def __repr__(self):
        return 'CaseColocalization: Cases ' + str(self)
        
    def valid_parameter_set(self, p_bounds=None, optimize=None, minimize=True, project=True, **kwargs):
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


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)