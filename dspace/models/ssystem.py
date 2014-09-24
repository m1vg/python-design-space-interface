''' Definition of the abstract model class.


'''

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.expressions import Expression

        
class SSystem(GMASystem):

    def __del__(self):
        if self._swigwrapper is not None:
            DSSSystemFree(self._swigwrapper)
    
    def set_swigwrapper(self, ssys_swigwrapper):
        self._swigwrapper = ssys_swigwrapper
        if self._swigwrapper is None:
            return
        Xd = VariablePool()
        Xd.set_swigwrapper(DSVariablePoolCopy(DSSSystemXd(ssys_swigwrapper)))
        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSSSystemXi(ssys_swigwrapper)))
        self._independent_variables = Xi
        self._dependent_variables = Xd
        
        
    def _parse_equations(self):
        auxiliary_variables = self.auxiliary_variables
        swigwrapper = DSSWIGSSystemParseWrapper(self.equations,
                                                len(self.equations),
                                                auxiliary_variables,
                                                len(auxiliary_variables)
                                                )
        self.set_swigwrapper(swigwrapper)
    
    @property
    def equations(self):
        sol = DSSSystemEquations(self._swigwrapper)
        solution = list()
        for i in xrange(0, DSSSystemNumberOfEquations(self._swigwrapper)):
            expr = DSExpressionAtIndexOfExpressionArray(sol, i)
            solution.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(sol)
        return Equations(solution, latex_symbols=self._latex)
    
    @property
    def solution(self):
        sol = DSSSystemSolution(self._swigwrapper)
        solution = list()
        for i in xrange(0, len(self.equations)):
            expr = DSExpressionAtIndexOfExpressionArray(sol, i)
            solution.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(sol)
        return Equations(solution, latex_symbols=self._latex)
    
    @property
    def solution_log(self):
        sol = DSSSystemLogarithmicSolution(self._swigwrapper)
        solution = list()
        for i in xrange(0, len(self.equations)):
            expr = DSExpressionAtIndexOfExpressionArray(sol, i)
            solution.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(sol)
        return Equations(solution, latex_symbols=self._latex)
    
    @property
    def dependent_variables(self):
        return self._dependent_variables.keys()
    
    def log_gain(self, dependent, independent):
        if dependent not in self.dependent_variables:
            raise NameError, str(dependent) + ' is not a dependent variable'
        if independent not in self.independent_variables:
            raise NameError, str(independent) + ' is not an independent variable'
        return DSSSystemLogarithmicGain(self._swigwrapper, dependent, independent)
    
    def remove_algebraic_constraints(self):
        return SSystem(self._equations,
                       name=self.name + '(ODE only)',
                       swigwrapper=DSSSystemByRemovingAlgebraicConstraints(self._swigwrapper))
        
    def steady_state(self, parameter_values, log_out=False):
        Xd = VariablePool()
        Xd.set_swigwrapper(DSSSystemXd(self._swigwrapper))
        steady_states = DSSSystemSteadyStateValues(self._swigwrapper, parameter_values._swigwrapper)
        var_names = Xd.keys()
        if log_out is False:
            steady_states = {var_names[i]:10**steady_states[i][0] for i in xrange(len(var_names))}
        else:
            steady_states = {var_names[i]:steady_states[i][0] for i in xrange(len(var_names))}
        Xd.set_swigwrapper(None)
        return steady_states
    
    def value_for_auxiliary_variables(self, Xdt0, Xi0, log_out=False):
        Xd_t = VariablePool()
        independent = VariablePool(names=self.independent_variables)
        for key in independent:
            value = float(Xi0[key])
            independent[key] = value
        for key in self.dependent_variables:
            if key in self.auxiliary_variables:
                continue
            value = float(Xdt0[key])
            Xd_t[key] = value      
        aux= DSSSystemAuxiliaryVariablesForSteadyState(self._swigwrapper,
                                                       Xd_t._swigwrapper,
                                                       independent._swigwrapper)
        j = 0
        if log_out == True:
            aux = [aux[i][0] for i in xrange(len(aux))]
        else:
            aux = [10**aux[i][0] for i in xrange(len(aux))]
        aux = {self.auxiliary_variables[i]:aux[i] for i in xrange(len(aux))}
        return aux

        
    def steady_state_flux(self, parameter_values, log_out=False):
        Xd = VariablePool()
        Xd.set_swigwrapper(DSSSystemXd(self._swigwrapper))
        flux = DSSSystemSteadyStateFlux(self._swigwrapper, parameter_values._swigwrapper)
        var_names = Xd.keys()
        if log_out is False:
            steady_states = {('V_' + var_names[i]):10**flux[i][0] for i in xrange(len(var_names))}
        else:
            steady_states = {('V_' + var_names[i]):flux[i][0] for i in xrange(len(var_names))}
        Xd.set_swigwrapper(None)
        return steady_states
    
    def steady_state_function(self, function, parameter_values):
        if isinstance(function, Expression):
            expr = function
        else:
            expr = Expression(function)
        p_vals = parameter_values.copy()
        p_vals.update(self.steady_state(parameter_values, log_out=False))
        p_vals.update(self.steady_state_flux(parameter_values, log_out=False))
        value = expr.eval_with_values(p_vals=p_vals)
        return value

    def positive_roots(self, parameter_values, show_marginal=True):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError, 'S-System must be reduced to ODE-only system'
        [positive_roots, has_marginal] = DSSSystemPositiveRootsSWIG(self._swigwrapper,
                                                                   parameter_values._swigwrapper)
        if has_marginal > 0:
            positive_roots = str(positive_roots) + '*'
        return positive_roots
    
    def routh_index(self, parameter_values):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError, 'S-System must be reduced to ODE only system'
        return DSSSystemRouthIndex(self._swigwrapper, parameter_values._swigwrapper)
        
    def routh_array(self, parameter_values):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError, 'S-System must be reduced to ODE only system'
        return DSSSystemRouthArray(self._swigwrapper, parameter_values._swigwrapper)
        
    
        
     