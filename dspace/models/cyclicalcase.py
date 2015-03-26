import itertools

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.models.case import Case, CaseIntersection
from dspace.models.ssystem import SSystem
from dspace.expressions import Expression

class CyclicalCase(Case):
    
    def __init__(self, equations, swigwrapper, name = None, latex_symbols=None, **kwargs):
        
        if name == None:
            name = 'Unnamed'
        super(Case, self).__init__(equations,
                                   name=name,
                                   latex_symbols=latex_symbols)
        setattr(self, '_ssystem', None)
        setattr(self, '_independent_variables', None)
        setattr(self, '_reduced_ssystem', None)
        setattr(self, '_freeData', False)
        self.set_swigwrapper(swigwrapper)
        
    def _cyclical_case(self, case, name):
        
        if isinstance(case, int) is False:
            raise TypeError, 'case must be indicated by its case number'
        sub=DSCyclicalCaseCyclicalSubcaseWithCaseNumber(self._swigwrapper, case)
        if sub is None:
            return None
        case = Case(self, DSCyclicalCaseSubcaseWithCaseNumber(self._swigwrapper, case), name)
        eq6=Equations(case.equations.system, case.auxiliary_variables)
        return CyclicalCase(eq6, sub, name = case.name)

    def __call__(self, index_or_iterable):
        if isinstance(index_or_iterable, (int, str)) is True:
            iterable = [index_or_iterable]
        else:
            iterable = index_or_iterable
        cases = list()
        for index in iterable:
            try:
                new_index = int(index)
                index = new_index
            except:
                pass
            if isinstance(index, int) is True:
                name = self.name + ': Subcase ' + str(index)
                case = self._cyclical_case(index, name)
                if case is None:
                    case = Case(self, DSCyclicalCaseSubcaseWithCaseNumber(self._swigwrapper, index), name=name)
                cases.append(case)
            elif isinstance(index, str) is True:
                indices = index.split('_')
                name = self.name + ': Subcase ' + indices[0]
                case = self._cyclical_case(int(indices[0]), name)
                if case is None:
                    raise ValueError, 'Subcase ' + str(indices[0]) + ' is not cyclical'
                indices.pop(0)
                last_index = int(indices[-1])
                indices.pop(-1)
                for i in indices:
                    name = name + ': Subcase ' + i
                    subindex = int(i)
                    case = case._cyclical_case(subindex, name)
                    if case is None:
                        raise ValueError, 'Subcase is not cyclical'
                name = name + ': Subcase ' + str(last_index)
                subcase = case(last_index)
                cases.append(subcase)
            else:
                raise TypeError, 'input argument must be a case number'
        if isinstance(index_or_iterable, (int, str)) is True:
            return cases[0]
        return cases
    
    def __del__(self):
        if self._freeData is True:
            DSCyclicalCaseFree(self._swigwrapper)
        return
    
    def __getstate__(self):
        odict = self.__dict__.copy()
        odict['_swigwrapper'] = DSSWIGDSCyclicalCaseEncodedBytes(self._swigwrapper)
        del odict['_ssystem']
        del odict['_independent_variables']
        del odict['_dependent_variables']
        return odict
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        encoded = state['_swigwrapper']
        self.set_swigwrapper(DSSWIGDSCyclicalCaseDecodeFromByteArray(encoded)) 
        self._freeData = True
        
    @property
    def equations(self):
        eqs = DSCyclicalCaseEquations(self._swigwrapper)
        equations = list()
        for i in xrange(0, DSCyclicalCaseNumberOfEquations(self._swigwrapper)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equations.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        return Equations(equations, latex_symbols=self._latex)
    
    @property
    def augmented_equations(self):
        eqs = DSDesignSpaceEquations(DSCyclicalCaseInternalDesignSpace(self._swigwrapper))
        equations = list()
        for i in xrange(0, DSCyclicalCaseNumberOfEquations(self._swigwrapper)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equations.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        return Equations(equations)
    
    @property
    def ssystem(self):
        return self._ssystem
    
    @property
    def independent_variables(self):
        return self._independent_variables.keys()
    
    @property
    def valid_parameter_set(self):
        subcases = self.valid_subcases()
        pvals_dict = dict()
        for i in subcases:
            pvals = self(i).valid_parameter_set
            pvals_dict[i] = pvals
        return pvals_dict
        
    def valid_interior_parameter_set(self, distance=50):
        subcases = self.valid_subcases()
        pvals_dict = dict()
        for i in subcases:
            pvals = self(i).valid_interior_parameter_set(distance=distance)
            pvals_dict[i] = pvals
        return pvals_dict

    @property
    def case_number(self):
        return DSCyclicalCaseIdentifier(self._swigwrapper)
    
    @property
    def signature(self):
        return DSCyclicalCaseSignatureToString(self._swigwrapper)
    
    @property
    def conditions(self):
        conditions = list()
        eqs_expr = DSCyclicalCaseConditions(self._swigwrapper)
        for i in xrange(0, DSCyclicalCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(conditions)
        
    @property
    def conditions_log(self):
        conditions = list()
        eqs_expr = DSCyclicalCaseLogarithmicConditions(self._swigwrapper)
        for i in xrange(0, DSCyclicalCaseNumberOfConditions(self._swigwrapper)):
            conditions.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        return Equations(conditions)
    
    @property
    def boundaries(self):
        return None
    
    @property
    def boundaries_log(self):
        return None
    
    @property
    def number_of_subcases(self):
        return DSCyclicalCaseNumberOfSubcases(self._swigwrapper)
        
    @property
    def original_case(self):
        case = Case(self, DSCaseCopy(DSCyclicalCaseOriginalCase(self._swigwrapper)), self.name + ' [original]')
        return case

    def set_swigwrapper(self, swigwrapper):
        self._swigwrapper = swigwrapper
        ## ds_swigwrapper = DSCyclicalCaseInternalDesignSpace(swigwrapper)
        
        Xd = VariablePool()
        Xd.set_swigwrapper(DSCyclicalCaseXd(swigwrapper))
        for i in VariablePool():
            if i not in self.dependent_variables:
                raise NameError, 'Dependent Variables are inconsistent'
        Xi = VariablePool()
        Xi.set_swigwrapper(DSCyclicalCaseXi(swigwrapper))
        self._independent_variables = Xi.copy()
        self._dependent_variables = Xd.copy()       
        Xd.set_swigwrapper(None)
        Xi.set_swigwrapper(None)
        eqs = list()
        eqs_expr = DSSSystemEquations(DSCyclicalCaseSSystem(self._swigwrapper))
        for i in xrange(0, DSCyclicalCaseNumberOfEquations(self._swigwrapper)):
            eqs.append(DSExpressionAsString(DSExpressionAtIndexOfExpressionArray(eqs_expr, i)))
            DSExpressionFree(DSExpressionAtIndexOfExpressionArray(eqs_expr, i))
        DSSecureFree(eqs_expr)
        self._ssystem = SSystem(self._equations,
                                name=self.name,
                                swigwrapper=DSSSystemCopy(DSCyclicalCaseSSystem(swigwrapper)))

    
    @property
    def dependent_variables(self):
        return self._dependent_variables.keys()
    
    @property
    def is_cyclical(self):
        return isinstance(self, CyclicalCase)
    
    def _valid_subcases_bounded(self, p_bounds):
        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        for i in lower:
            lower[i] = 1E-20
            upper[i] = 1E20
        for (key,value) in p_bounds.iteritems():
            if key not in lower:
                continue
            try:
                min_value,max_value = value
            except TypeError:
                min_value = value
                max_value = value
            if min_value > max_value:
                raise ValueError, 'parameter slice bounds are inverted: min is larger than max'
            lower[key] = min_value
            upper[key] = max_value
        valid_cases = DSCyclicalCaseCalculateAllValidSubcasesForSlice(self._swigwrapper,
                                                                      lower._swigwrapper,
                                                                      upper._swigwrapper)
        number_of_cases = DSDictionaryCount(valid_cases)
        cases = list()
        keys = [DSDictionaryKeyAtIndex(valid_cases, i) for i in xrange(0, number_of_cases)]
        for key in keys:
            case_swigwrapper = DSSWIGVoidAsCase(DSDictionaryValueForName(valid_cases, key))
            cases.append(key.split('_')[1])
            DSCaseFree(case_swigwrapper)
        DSDictionaryFree(valid_cases)
        cases.sort()
        return cases
       
    def valid_subcases(self, p_bounds=None):
        if p_bounds is not None:
            return self._valid_subcases_bounded(p_bounds)
        valid_cases = DSCyclicalCaseCalculateAllValidSubcases(self._swigwrapper)
        number_of_cases = DSDictionaryCount(valid_cases)
        cases = list()
        keys = [DSDictionaryKeyAtIndex(valid_cases, i) for i in xrange(0, number_of_cases)]
        for key in keys:
            case_swigwrapper = DSSWIGVoidAsCase(DSDictionaryValueForName(valid_cases, key))
            cases.append(key.split('_')[1])
            DSCaseFree(case_swigwrapper)
        DSDictionaryFree(valid_cases)
        cases.sort()
        return cases
    
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
        dictionary=DSCyclicalCaseVerticesFor2DSlice(self._swigwrapper, 
                                                    lower._swigwrapper,
                                                    upper._swigwrapper,
                                                    x_variable,
                                                    y_variable)
        number_of_cases = DSDictionaryCount(dictionary)
        all_vertices = dict()
        keys = [DSDictionaryKeyAtIndex(dictionary, i) for i in xrange(0, number_of_cases)]
        for key in keys:
            log_vertices = DSSWIGVoidAsVertices(DSDictionaryValueForName(dictionary, key))
            if log_out is True:
                vertices = log_vertices
            else:
                vertices = list()
                for vertex in log_vertices:
                    vertices.append([10**coordinate for coordinate in vertex])
            all_vertices[key.split('_')[1]] = vertices
        DSDictionaryFree(dictionary)
        return all_vertices
    
    def steady_state(self, parameter_values, log_out=False):
        Xd = VariablePool(names=self.dependent_variables)
        p_vals = VariablePool(names=self.independent_variables)
        for i in p_vals:
            p_vals[i] = parameter_values[i]
        cases = self.valid_subcases(p_bounds=p_vals)
        steady_states = dict()
        for i in cases:
            case = self(i)
            ssys = case.ssystem
            ss = DSSSystemSteadyStateValues(ssys._swigwrapper, parameter_values._swigwrapper)
            var_names = Xd.keys()
            if log_out is False:
                steady_states[str(i)] = {var_names[j]:10**ss[j][0] for j in xrange(len(var_names))}
            else:
                steady_states[str(i)] = {var_names[j]:ss[j][0] for j in xrange(len(var_names))}
        return steady_states
        
    def steady_state_flux(self, parameter_values, log_out=False):
        steady_states = self.steady_state(parameter_values)
        Xd = VariablePool(names=self.dependent_variables)
        fluxes = dict()
        case_swigwrapper = DSCyclicalCaseOriginalCase(self._swigwrapper)
        for i in steady_states:
            Xd.update(steady_states[i])
            flux = DSSSystemSteadyStateFluxForDependentVariables(DSCaseSSystem(case_swigwrapper),
                                                                               Xd._swigwrapper,
                                                                               parameter_values._swigwrapper)
            var_names = Xd.keys()
            if log_out is False:
                flux = {('V_' + var_names[j]):10**flux[j][0] for j in xrange(len(var_names))}
            else:
                flux = {('V_' + var_names[j]):flux[j][0] for j in xrange(len(var_names))}
            fluxes[i] = flux
        return fluxes
    
    def steady_state_function(self, function, parameter_values):
                
        if isinstance(function, Expression):
            expr = function
        else:
            expr = Expression(function)
        p_values = VariablePool(names=self.independent_variables)
        for key in p_values:
            p_values[key] = parameter_values[key]
        ss = self.steady_state(p_values, log_out=False)
        flux = self.steady_state_flux(p_values, log_out=False)
        values = dict()
        for subcase in ss:
            p_values.update(ss[subcase])
            p_values.update(flux[subcase])
            value = expr.eval_with_values(p_vals=p_values)
            values[subcase] = value
        return values
        
    def positive_roots(self, parameter_values):
        steady_states = self.steady_state(parameter_values)
        fluxes = self.steady_state_flux(parameter_values)
        roots = dict()
        case_swigwrapper = DSCyclicalCaseOriginalCase(self._swigwrapper)
        if self._reduced_ssystem is None:
            self._reduced_ssystem = DSSSystemByRemovingAlgebraicConstraints(DSCaseSSystem(DSCyclicalCaseOriginalCase(self._swigwrapper)))
        for i in steady_states:
            case = self(i)
            ssys = case.ssystem.remove_algebraic_constraints()
            #Xd = VariablePool(names=ssys.dependent_variables)
            #flux = VariablePool(names=['V_' + j for j in ssys.dependent_variables])
            #Xd.update({j:steady_states[i][j] for j in ssys.dependent_variables})
            #flux.update({('V_'+j):fluxes[i]['V_'+j] for j in ssys.dependent_variables})
            #flux.update(ssys.steady_state_flux(parameter_values))
            #Xd.update(ssys.steady_state(parameter_values))
            #roots[i] = DSSSystemPositiveRootsForSteadyStateAndFlux(ssys._swigwrapper,
            #                                                       Xd._swigwrapper,
            #                                                       parameter_values._swigwrapper,
            #                                                       flux._swigwrapper)
            roots[i] = ssys.positive_roots(parameter_values)
        return roots
    
    def _parse_equations(self):
        return    
        
