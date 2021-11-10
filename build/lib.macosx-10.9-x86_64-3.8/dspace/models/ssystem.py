''' Definition of the abstract model class.


'''

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.expressions import Expression

import numpy as np
import time
        
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
        Xd_t = VariablePool()
        Xd_t.set_swigwrapper(DSVariablePoolCopy(DSSSystemXd_t(ssys_swigwrapper)))
        if DSSSystemXd_a_c(ssys_swigwrapper) is not None:
            Xd_a_c = VariablePool()
            Xd_a_c.set_swigwrapper(DSVariablePoolCopy(DSSSystemXd_a_c(ssys_swigwrapper)))
        else:
            Xd_a_c = None
        self._independent_variables = Xi
        self._dependent_variables = Xd
        self._dependent_variables_no_algebraic = Xd_t
        self._conserved_variables = Xd_a_c
        
        
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
        for i in range(0, DSSSystemNumberOfEquations(self._swigwrapper)):
            expr = DSExpressionAtIndexOfExpressionArray(sol, i)
            solution.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(sol)
        return Equations(solution, latex_symbols=self._latex)
    
    @property
    def m(self):
        return DSSSystemM(self._swigwrapper)
    
    @property
    def Ad(self):
        return DSSSystemAd(self._swigwrapper)
    
    @property
    def Ai(self):
        return DSSSystemAi(self._swigwrapper)
    
    @property
    def Gd(self):
        return DSSSystemGd(self._swigwrapper)

    @property
    def Hd(self):
        return DSSSystemHd(self._swigwrapper)

    @property
    def alpha(self):
        return DSSSystemAlpha(self._swigwrapper)

    @property
    def beta(self):
        return DSSSystemBeta(self._swigwrapper)

    @property
    def Gi(self):
        return DSSSystemGi(self._swigwrapper)

    @property
    def Hi(self):
        return DSSSystemHi(self._swigwrapper)

    @property
    def solution(self):
        if DSSSystemHasSolution(self._swigwrapper) is False:
            if DSSSystemIsUnstable(self._swigwrapper) is True:
                sol = DSuSSystemSolution(self._swigwrapper)
                solution = list()
                dependent_var = DSVariablePoolNumberOfVariables(DSSSystemXd_t(self._swigwrapper))
                for i in range(0, dependent_var):
                    expr = DSExpressionAtIndexOfExpressionArray(sol, i)
                    solution.append(DSExpressionAsString(expr))
                    DSExpressionFree(expr)
                DSSecureFree(sol)
                return Equations(solution, latex_symbols=self._latex)
            else:
                return None
        else:
            sol = DSSSystemSolution(self._swigwrapper)
            solution = list()
            for i in range(0, len(self.equations)):
                expr = DSExpressionAtIndexOfExpressionArray(sol, i)
                solution.append(DSExpressionAsString(expr))
                DSExpressionFree(expr)
            DSSecureFree(sol)
            return Equations(solution, latex_symbols=self._latex)
    
    @property
    def solution_log(self):
        if DSSSystemHasSolution(self._swigwrapper) is False:
            if DSSSystemIsUnstable(self._swigwrapper) is True:
                sol = DSuSSystemLogarithmicSolution(self._swigwrapper)
                solution = list()
                dependent_var = DSVariablePoolNumberOfVariables(DSSSystemXd_t(self._swigwrapper))
                for i in range(0, dependent_var):
                    expr = DSExpressionAtIndexOfExpressionArray(sol, i)
                    solution.append(DSExpressionAsString(expr))
                    DSExpressionFree(expr)
                DSSecureFree(sol)
                return Equations(solution, latex_symbols=self._latex)
            else:
                return None
        else:
            sol = DSSSystemLogarithmicSolution(self._swigwrapper)
            solution = list()
            for i in range(0, len(self.equations)):
                expr = DSExpressionAtIndexOfExpressionArray(sol, i)
                solution.append(DSExpressionAsString(expr))
                DSExpressionFree(expr)
            DSSecureFree(sol)
            return Equations(solution, latex_symbols=self._latex)
    
    @property
    def dependent_variables(self):
        return self._dependent_variables.keys()
    @property
    def dependent_variables_no_algebraic(self):
        return self._dependent_variables_no_algebraic.keys()
    @property
    def conserved_variables(self):
        if self._conserved_variables is not None:
            return self._conserved_variables.keys()
        else:
            return []
    @property
    def should_adjust_codominant_stoichiometry(self):
        return DSSSystemAdjustCodominantStoichiometry(self._swigwrapper)
    
    def log_gain(self, dependent, independent):
        if dependent not in self.dependent_variables:
            raise NameError(str(dependent) + ' is not a dependent variable')
        if independent not in self.independent_variables:
            raise NameError(str(independent) + ' is not an independent variable')
        if DSSSystemIsUnstable(self._swigwrapper) is False:
            return DSSSystemLogarithmicGain(self._swigwrapper, dependent, independent)
        else:
            return DSuSSystemLogarithmicGain(self._swigwrapper, dependent, independent)
    
    def remove_algebraic_constraints(self):
        return SSystem(self._equations,
                       name=self.name + '(ODE only)',
                       swigwrapper=DSSSystemByRemovingAlgebraicConstraints(self._swigwrapper),
                       latex_symbols=self._latex)
        
    def steady_state(self, parameter_values, log_out=False):
        if DSSSystemIsUnstable(self._swigwrapper) is True:
            return self.steady_state_blowing(parameter_values, log_out)
        if DSSSystemHasSolution(self._swigwrapper) is False:
            return None
        if DSSSystemAdjustCodominantStoichiometry(self._swigwrapper) is True:
            DSSSystemAdjustStoichiometryOfCodominantCase(self._swigwrapper)
        Xd = VariablePool()
        Xd.set_swigwrapper(DSSSystemXd(self._swigwrapper))
        steady_states = DSSSystemSteadyStateValues(self._swigwrapper, parameter_values._swigwrapper)
        var_names = Xd.keys()
        if log_out is False:
            steady_states = {var_names[i]:10**steady_states[i][0] for i in range(len(var_names))}
        else:
            steady_states = {var_names[i]:steady_states[i][0] for i in range(len(var_names))}
        Xd.set_swigwrapper(None)
        return steady_states

    def steady_state_blowing(self, parameter_values, log_out=False):
        Xd_t = VariablePool()
        Xd_t.set_swigwrapper(DSSSystemXd_t(self._swigwrapper))
        steady_states = DSuSSystemSteadyStateValues(self._swigwrapper, parameter_values._swigwrapper)
        var_names = Xd_t.keys()
        if log_out is False:
            steady_states = {var_names[i]:10**steady_states[i][0] for i in range(len(var_names))}
        else:
            steady_states = {var_names[i]:steady_states[i][0] for i in range(len(var_names))}
        Xd_t.set_swigwrapper(None)

        ## This section is used to add values for Xd_a_c to the dictionary steady_states.
        if DSSSystemIsConserved(self._swigwrapper) is True:
            Xd_a_c = VariablePool()
            Xd_a_c.set_swigwrapper(DSSSystemXd_a_c(self._swigwrapper))
            steady_states_conserved = DSuSSystemSteadyStateValuesForConservedVariables(self._swigwrapper, parameter_values._swigwrapper)
            if steady_states_conserved is None:
                return steady_states
            var_names = Xd_a_c.keys()
            if log_out is False:
                steady_states.update({var_names[i]: 10 ** steady_states_conserved[i][0] for i in range(len(var_names))})
            else:
                steady_states.update({var_names[i]: steady_states_conserved[i][0] for i in range(len(var_names))})
            Xd_a_c.set_swigwrapper(None)
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
            aux = [aux[i][0] for i in range(len(aux))]
        else:
            aux = [10**aux[i][0] for i in range(len(aux))]
        aux = {self.auxiliary_variables[i]:aux[i] for i in range(len(aux))}
        return aux

        
    def steady_state_flux(self, parameter_values, log_out=False):
        Xd = VariablePool()
        Xd.set_swigwrapper(DSSSystemXd_t(self._swigwrapper))
        if DSSSystemAdjustCodominantStoichiometry(self._swigwrapper) is True:
            DSSSystemAdjustStoichiometryOfCodominantCase(self._swigwrapper)
        if DSSSystemIsUnstable(self._swigwrapper) is True:
            flux = DSuSSystemSteadyStateFlux(self._swigwrapper, parameter_values._swigwrapper)
        else:
            flux = DSSSystemSteadyStateFlux(self._swigwrapper, parameter_values._swigwrapper)
        var_names = Xd.keys()
        if log_out is False:
            steady_states = {('V_' + var_names[i]):10**flux[i][0] for i in range(len(var_names))}
        else:
            steady_states = {('V_' + var_names[i]):flux[i][0] for i in range(len(var_names))}
        Xd.set_swigwrapper(None)

        # ## This section is used to add values for Xd_a_c to the dictionary of fluxes.
        # if DSSSystemIsConserved(self._swigwrapper) is True:
        #     Xd_a_c = VariablePool()
        #     Xd_a_c.set_swigwrapper(DSSSystemXd_a_c(self._swigwrapper))
        #     flux_conserved = DSSSystemSteadyStateFluxForConservedVariables(self._swigwrapper, parameter_values._swigwrapper)
        #     var_names = Xd_a_c.keys()
        #     if log_out is False:
        #         steady_states.update({var_names[i]: 10 ** flux_conserved[i][0] for i in xrange(len(var_names))})
        #     else:
        #         steady_states.update({var_names[i]: flux_conserved[i][0] for i in xrange(len(var_names))})
        #     Xd_a_c.set_swigwrapper(None)

        return steady_states
    
    def steady_state_function(self, function, parameter_values):
        if isinstance(function, Expression):
            expr = function
        else:
            expr = Expression(function)
        if DSSSystemAdjustCodominantStoichiometry(self._swigwrapper) is True:
            DSSSystemAdjustStoichiometryOfCodominantCase(self._swigwrapper)
        p_vals = parameter_values.copy()
        p_vals.update(self.steady_state(parameter_values, log_out=False))
        p_vals.update(self.steady_state_flux(parameter_values, log_out=False))
        for i in self.dependent_variables_no_algebraic:
            for j in self.independent_variables:
                p_vals['$L_'+i+'_'+j] = self.log_gain(i,j)
        value = expr.eval_with_values(p_vals=p_vals)
        return value

    def positive_roots(self, parameter_values, show_marginal=True):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError('S-System must be reduced to ODE-only system')
        [positive_roots, has_marginal] = DSSSystemPositiveRootsSWIG(self._swigwrapper,
                                                                    parameter_values._swigwrapper)
        if has_marginal > 0:
            positive_roots = str(positive_roots) + '*'
        return positive_roots
        
    def eigenvalues(self, parameter_values):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError('S-System must be reduced to ODE-only system')
        ss = self.steady_state(parameter_values)
        fluxes = self.steady_state_flux(parameter_values)
        turnover = np.array(zip([fluxes['V_'+i]/ss[i] for i in self.dependent_variables]))
        FAd = turnover * np.array(self.Ad)
        eigenvalues = np.linalg.eig(FAd)[0]
        return eigenvalues

    def positive_roots_numpy(self, parameter_values):
        eigen_values = self.eigenvalues(parameter_values)
        num_pos = len([x for x in eigen_values.real if x > 0])
        return num_pos

    def has_complex_conjugates(self, parameter_values):
        eigen_values = self.eigenvalues(parameter_values)
        has_conjugates = False

        # count_zero = 0
        # for idx, x in enumerate(eigen_values.imag):
        #     if x == 0:
        #         count_zero += 1
        # if count_zero == len(eigen_values):
        #     return has_conjugates

        for val1 in eigen_values.imag:
            if has_conjugates is True:
                break
            for val2 in eigen_values.imag:
                if val1 == - val2:
                    if val1 != 0 and val2 != 0:
                        has_conjugates = True
                        break

        return has_conjugates
    
    def routh_index(self, parameter_values):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError('S-System must be reduced to ODE only system')
        return DSSSystemRouthIndex(self._swigwrapper, parameter_values._swigwrapper)
        
    def routh_array(self, parameter_values):
        
        if DSVariablePoolNumberOfVariables(DSSSystemXd_a(self._swigwrapper)) > 0:
            raise TypeError('S-System must be reduced to ODE only system')
        return DSSSystemRouthArraySWIG(self._swigwrapper, parameter_values._swigwrapper)
        
    
        
     