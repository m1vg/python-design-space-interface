import itertools
import sys
import dspace
from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.models.ssystem import SSystem
from dspace.models.case import Case, CaseIntersection, CaseColocalization
from dspace.models.cyclicalcase import CyclicalCase
from dspace.expressions import Expression
from dspace.models.designspace import DesignSpace
import numpy as np
from scipy.integrate import odeint
import time
import numpy as np
from functools import cmp_to_key
import difflib
from IPython.utils import io
with io.capture_output() as capture:
    from assimulo.solvers import IDA
    from assimulo.solvers import Radau5DAE
    from assimulo.problem import Implicit_Problem

class FullSystem(object):

    def __init__(self, ds):

        setattr(self, 'ds', ds)

        self.coef_dic = {}
        self.conservations_dictionary = None

        # 1. Get the auxiliary, independent and dependent variables
        auxiliary = ds.auxiliary_variables
        dependent = ds.dependent_variables
        dependent_no_artificial = ds.dependent_variables

        # 2. Remove auxiliary variables from dependent variables
        for aux in auxiliary:
            if aux in dependent:
                dependent.remove(aux)

        # 3. Generate conserved, artificial variables Xc1, ... Xcn
        if ds.number_conservations != 0:
            artificial_var = ['Xc' + str(n + 1) for n in range(int(ds.number_conservations))]

            # Remove artificial variables from dependent variables
            for aux in artificial_var:
                if aux in dependent_no_artificial:
                    dependent_no_artificial.remove(aux)

        setattr(self, 'dependent_no_aux', dependent)
        self.generate_ode_equations()
        dependent_no_aux_no_conserved = [variable for variable in dependent if variable not in self.coef_dic.keys()]
        setattr(self, 'dependent_no_aux_no_conserved', dependent_no_aux_no_conserved)
        self.dependent_no_artificial = dependent_no_artificial

    def generate_ode_equations(self):

        # construct a list  of len(dependent) first elements of ds.equations.
        depend_expression = self.ds.equations[0:len(self.dependent_no_aux)]
        number_conservations = self.ds.number_conservations
        n = len(self.dependent_no_aux)
        aux_expression = self.ds.equations[n:
                                           n + len(self.ds.auxiliary_variables) - number_conservations]

        # use the subst method of the expression object to substitute auxiliary expressions into dependent
        # expressions. first, generate a dictionary out of aux_expression list.
        aux_dic = dict((str(element.lhs), element.rhs) for element in aux_expression)

        # use the subst method on the dictionary itself to take care of nested auxiliary variables present in metabolic models.
        for element in aux_dic:
            aux_dic[element] = str(aux_dic[element].subst(aux_dic))
        self.aux_dic = aux_dic

        # substitute and create list.
        if len(self.ds.auxiliary_variables) - number_conservations != 0:
            depend_expression = [element.subst(aux_dic) for element in depend_expression]

        # check if there are conservation constraints present. If so, modify the array for depend_expression
        if number_conservations != 0:
            self.conservations = self.ds.equations[-number_conservations:]
            depend_expression = self.integrate_conserved_relationships(depend_expression)
            self.dae_equations = depend_expression + self.conservations
        self.ode_equations = depend_expression

    def integrate_conserved_relationships(self, depend_expression):

        delete_list = []
        modified_depend_expression = []
        for conservation in self.conservations:

            # first get the list of dependent variables involved in the conservation
            involved_conserved_variables = [var for var in conservation.variables if var in self.dependent_no_aux]
            unique_involved_conserved_variables = [var for var in involved_conserved_variables if
                                                   var not in delete_list]
            if len(unique_involved_conserved_variables) != 0:
                var_to_delete = unique_involved_conserved_variables[0]
                delete_list.append(var_to_delete)
            else:
                print("Please Check Syntax of Conservation Constraints")

        # Now loop over depend_expression. Delete if left hand side is equal to var_to_delete
        for equation in depend_expression:
            if str(equation.lhs).replace('.', '') not in delete_list:
                modified_depend_expression.append(equation)
        self.process_conservation_equations(delete_list)
        return modified_depend_expression

    def process_conservation_equations(self, variables):
        # variables is a list of string with variables matching self.conservations.
        self.conservations_dictionary = {}
        self.conservations_equal_zero_dictionary = {}
        self.conservations_variable_list = variables

        if len(variables) != len(self.conservations):
            print("Warning! Conservation relationships won't be correctly parsed. Don't trust results.")

        for indx, eq in enumerate(self.conservations):
            parts = str(eq).split("-")

            # 1. We first verify that the variables[index] is present in the respective conservation and extract the
            # expression in which the variables[index] is present. This is stored in element_to_replace
            element_to_replace = ''
            for element in parts:
                if variables[indx] in dspace.Expression(element).variables:
                    parts_to_replace = variables[indx]
                    element_to_replace = element

            # print("The variable[index] is ", variables[indx])
            # print("The element_to_replace is ", element_to_replace)

            # 2. We check if element_to_replace is different from variables[index]. If so, we need to extract the multipliers
            if element_to_replace == '':
                print("Warning! Conservation relationships won't be correctly parsed. Don't trust results.")
            else:
                if element_to_replace == variables[indx]:
                    multiplier = ''
                else:
                    multiplier = self.get_multipliers(variables[indx], element_to_replace)

            # print("The multipliers are: ", multiplier)

            # 3. Loop over elements of parts, exclude element_to_replace and add multiplier.
            parts_modified = []
            for element in parts:
                if element != element_to_replace:
                        parts_modified.append(element + multiplier)

            # print("parts_modified is: ", parts_modified)
            # print("parts_to_replace is: ", parts_to_replace)

            # 3. We loop over parts_modified and construct the subset parts_filtered, which excludes the element_to_replace
            parts_filtered =[]
            for element in parts_modified:
                parts_filtered.append(element.replace("0", parts_to_replace))

            # print("parts_filtered is ", parts_filtered)

            if len(parts_to_replace) == 0:
                print("the len of parts_to_replace is ", len(parts_to_replace))
                print("Warning, parsing of algebraic constraints is wrong! Do not trust results.")
                print("parts_to_replace  is ", parts_to_replace)
                print("variables are ", variables)

            new_eq = "-".join(parts_filtered)
            # print("The new equation is: ", new_eq)

            self.conservations[indx] = Expression(new_eq)
            self.conservations_dictionary.update({variables[indx]: Expression(new_eq)})
            self.conservations_equal_zero_dictionary.update({variables[indx]: eq})

            # Now we construct a factor dicitonary to take care of coefficients in the LHS
            # Note that the dictionary self.coef_dic is not required anymore.

            coef_list = parts_to_replace[0].split(variables[indx])
            try:
                coef = float(coef_list[0])
            except:
                coef = 1
            self.coef_dic.update({variables[indx]: coef})

    def get_multipliers(self, s, s_mult):

        output_list = [li for li in difflib.ndiff(s, s_mult) if li[0] != ' ']
        processed_output_list = []
        for element in output_list:
            processed_output_list.append(element.replace('+ ', ''))
        multiplier = "".join(processed_output_list)

        ## delete a * sign at the beginning
        if multiplier[0] == '*':
            multiplier = multiplier[1:]

        # delete * sign in the last position.
        if multiplier[-1] == '*':
            multiplier = multiplier[:-1]

        # replace double ** by single
        multiplier = multiplier.replace('**', '*')

        # split multiplier in multiple elements
        multiplier = multiplier.split('*')

        # loop over each element and replace ^ by ^- and add ^-1 if ^ not present:
        multiplier_mod = []
        for l in multiplier:
            if l.find('^') != -1:
                # replace ^ sign by ^-
                multiplier_mod.append(l.replace('^', '^-'))
            else:
                multiplier_mod.append(l + '^-1')

        multiplier = '*'.join(multiplier_mod)

        # delete * sign in the last position.
        if multiplier[-1] == '*':
            multiplier = multiplier[:-1]
        # add * sign in the first position
        if multiplier[0] != '*':
            multiplier = '*' + multiplier
        return multiplier

    def time_course_case(self, pvals, to, tf, nt, initial_concentrations=None, default_xo=0.001):

        # create parameters pool & make it available for all methods.
        parameters = dspace.VariablePool(names=self.dependent_no_aux + self.ds.independent_variables)
        parameters.update(pvals)
        self.parameters = parameters
        t = np.linspace(to, tf, nt)

        if initial_concentrations is None:
            initial_concentrations = dict((i, default_xo) for i in self.dependent_no_aux)

        # set initial conditions.
        yinit = np.array([])
        yinit_conserved = np.array([])

        for element in self.ode_equations:
            yinit = np.append(yinit, initial_concentrations[str(element.lhs).replace('.', '')])

        if self.ds.number_conservations != 0:
            for element in self.conservations_variable_list:
                yinit_conserved = np.append(yinit_conserved, initial_concentrations[element])

        results_dic = self.solve_ode(yinit, yinit_conserved, t)
        self.complement_results_dic_with_aux_variables(results_dic, len(t), pvals)

        return t, results_dic

    def titrations_case(self, pvals, to, tf, nt, variable, initial_concentrations=None, default_xo=0.001,
                        titration_range_min=10**-3, titration_range_max=10**3, variable_steps=100,
                        solver='ODEINT'):

        # create parameters pool & make it available for all methods.
        parameters = dspace.VariablePool(names=self.dependent_no_aux + self.ds.independent_variables)
        parameters.update(pvals)
        self.parameters = parameters
        t = np.linspace(to, tf, nt)

        if initial_concentrations is None:
            initial_concentrations = dict((i, default_xo) for i in self.dependent_no_aux)

        yinit = np.array([])
        yinit_conserved = np.array([])

        for element in self.ode_equations:
            yinit = np.append(yinit, initial_concentrations[str(element.lhs).replace('.', '')])

        if self.ds.number_conservations != 0:
            for element in self.conservations_variable_list:
                yinit_conserved = np.append(yinit_conserved, initial_concentrations[element])

        # Generate vector containing values for the selected variable on the x-axis.
        target_x_parameter = np.linspace(np.log10(titration_range_min),
                                         np.log10(titration_range_max),
                                         num=variable_steps)

        target_x_parameter = 10 ** target_x_parameter

        # Actualize parameters
        self.parameters.update({variable: target_x_parameter[0]})

        # Generate first steady state point
        y = self.solve_ode(yinit, yinit_conserved, t)

        # parse solution.
        # create a dictionary with data.
        results_dic = {}
        for counter, element in enumerate(self.ode_equations):
            v = str(element.lhs).replace('.', '')
            results_dic.update({v: [y[v][-1]]})

        if solver != 'ODEINT' and self.ds.number_conservations != 0:
            for counter2, element in enumerate(self.conservations_variable_list):
                # results_dic.update({element: [y[-1, counter + counter2 + 1]]})
                results_dic.update({element: [y[element][-1]]})

        # loop over other elements of target_x_parameter
        for w in target_x_parameter[1:]:
            self.parameters.update({variable: w})
            # initial concentrations correspond to last steady state point in results_dic.
            # Remember that the order of the dictionary does not correspond to the order with wich the elements were stored.
            yinit = np.array([])
            yinit_conserved = np.array([])

            for element in self.ode_equations:
                yinit = np.append(yinit, results_dic[str(element.lhs).replace('.', '')][-1])

            if solver != 'ODEINT' and self.ds.number_conservations != 0:
                for element in self.conservations_variable_list:
                    yinit_conserved = np.append(yinit_conserved, results_dic[element][-1])

            # Generate next steady state point
            y = self.solve_ode(yinit, yinit_conserved, t)

            # parse and append
            for counter, element in enumerate(self.ode_equations):
                # results_dic[str(element.lhs).replace('.', '')].append(y[-1, counter])
                v = str(element.lhs).replace('.', '')
                results_dic[v].append(y[v][-1])
            if solver != 'ODEINT' and self.ds.number_conservations != 0:
                for counter2, element in enumerate(self.conservations_variable_list):
                    # results_dic[element].append(y[-1, counter + counter2 + 1])
                    results_dic[element].append(y[element][-1])

        if self.conservations_dictionary is not None and solver == 'ODEINT':
            self.complement_results_dic(results_dic)

        # complement dictionary with auxiliary variables
        self.complement_results_dic_with_aux_variables(results_dic, len(target_x_parameter), self.parameters)

        return target_x_parameter, results_dic

    def complement_results_dic_with_aux_variables(self, results_dic, n, pvals):

        # This function should add vectors for auxiliary variables (algebraic constraints) to the results_dictionary.
        # Expressions defining these variables are defined in self.aux_dic
        aux_dic = self.aux_dic
        # 1. update results_dic with keys of aux_dic:
        keys_aux_dic = aux_dic.keys()
        keys_dep_var = results_dic.keys()
        results_dic.update({key: [] for key in keys_aux_dic})
        p_values = pvals.copy()

        # 2. loop over vectors, update dictionary, compute value of auxiliary variables and update corresponding vectors
        for ndx in range(n):
            # update parameter dictionary with values for dependent variables
            p_values.update({key: results_dic[key][ndx] for key in keys_dep_var})
            # compute value of auxiliary variables and update corresponding vectors
            for element in keys_aux_dic:
                results_dic[element].append(Expression(aux_dic[element]).eval_with_values(p_vals=p_values))
        # 3. Transform list into np.arrays for auxiliary variables. This is done for data consistency.
        for element in keys_aux_dic:
            results_dic[element] = np.array(results_dic[element])

    def complement_results_dic(self, results_dic):

        conservation_dic = self.conservations_dictionary
        nr_rows = len(results_dic[results_dic.keys()[0]])
        nr_colums = len(conservation_dic.keys())
        new_results = np.zeros([nr_rows, nr_colums])
        parameters_dic = self.parameters.copy()

        for row, _ in enumerate(results_dic[results_dic.keys()[0]]):
            #update dictionary of parameters
            sub_dic = dict((var, results_dic[var][row]) for var in results_dic.keys())
            parameters_dic.update(sub_dic)
            for col, key in enumerate(conservation_dic.keys()):
                expression = conservation_dic[key]
                new_results[row,col] = expression.rhs.eval_with_values(p_vals=parameters_dic) / self.coef_dic[key]

        # now update results_dic with rows from new_results matrix

        for col, key in enumerate(conservation_dic.keys()):
            results_dic.update({key: new_results[:,col]})

    def solve_ode(self, yinit, yinit_conserved, time, solver='ODEINT'):

        ode = self.ode_equations
        n_conservations = self.ds.number_conservations
        if n_conservations != 0:
            conservations = self.conservations_equal_zero_dictionary
            conservations_list = self.conservations_variable_list

        # Standard ODE, no conservations.
        if n_conservations == 0:
            # rtol = 1e-10  #defaults are 1.49012e-8
            # atol = 1e-10  #defaults are 1.49012e-8
            y = odeint(self.construct_fun, yinit, time)
            results_dic = {}
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                results_dic.update({var: y[:, counter]})
            return results_dic

        # Recommended option to solve DAE systems.
        elif n_conservations != 0 and solver != 'ODEINT':

            # 0. Expand y0 to include algebraic constraints
            y0 = np.append(yinit, yinit_conserved)

            # 1. calculate yd0 from yd
            yd0 = self.calculate_yd0_from_y0(yinit, yinit_conserved)

            # 2. create implicit problem
            imp_mod = Implicit_Problem(self.construct_fun_DAE, y0, yd0)

            # 3. set the algebraic components
            v = np.zeros(len(ode) + n_conservations) + 1
            v[-n_conservations:] = 0
            imp_mod.algvar = v

            if solver == 'IDA':
                # Create an Assimulo implicit solver (IDA)
                imp_sim = IDA(imp_mod)  # Create a IDA solver

                # Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
                imp_sim.make_consistent('IDA_YA_YDP_INIT')

            elif solver == 'RADAU5':
                imp_sim = Radau5DAE(imp_mod)

            # Sets the paramters
            imp_sim.atol = 1e-6  # Default 1e-6
            imp_sim.rtol = 1e-6  # Default 1e-6
            imp_sim.suppress_alg = False  # Suppress the algebraic variables on the error test
            imp_sim.verbosity = 50 # Default 50, scream = 10
            imp_sim.maxsteps = 10000  # default 10000

            # Simulate
            t, y, yd = imp_sim.simulate(time[-1], len(time)-1)
            if self.type.value == 'Titrations':
                return y

            # pack results in a dictionary. First get vectors for ode_variables
            results_dic = {}
            for counter, element in enumerate(self.ode_equations):
                var = str(element.lhs).replace('.', '')
                results_dic.update({var: y[:, counter]})
            # now get vectors for algebraic constraints
            for counter2, element in enumerate(conservations_list):
                n = counter + counter2 + 1
                results_dic.update({element: y[:, n]})
            # if len(self.advanced_time_courses.ss_functions.members) > 0:
            #     self.calculate_ss_function(results_dic)
            return results_dic

    def construct_fun(self, y, t):

        # initialize the output array
        y_out = np.array([])

        # Parameter values of this object are actualized to reflect parameter values of the main toolbox.
        for counter, element in enumerate(self.ode_equations):
            self.parameters[str(element.lhs).replace('.', '')] = y[counter]

        # if conserved relationships are present, update the dictionary self.parameters with the respective values.
        if self.conservations_dictionary is not None:
            for conserved_var, expression in self.conservations_dictionary.items():
                self.parameters[conserved_var] = expression.rhs.eval_with_values(p_vals=self.parameters) / self.coef_dic[conserved_var]

        for element in self.ode_equations:
            y_out = np.append(y_out, element.rhs.eval_with_values(p_vals=self.parameters))

        return y_out

    def calculate_yd0_from_y0(self, yinit, yinit_conserved):

        yd0 = []

        # create dictionary of parameters and update it by taking values from toolbox
        par = self.parameters

        # update dictionary of parameter values to consider initial conditions for Xd_t pools
        for indx, element in enumerate(self.ode_equations):
            par[str(element.lhs).replace('.', '')] = yinit[indx]

        # update dictionary of parameter values to consider intial conditions for Xd_a_c pools
        for indx, element in enumerate(self.conservations_variable_list):
            par[element] = yinit_conserved[indx]

        # Now start populating yd0. Start evaluating the rhs of self.ode_equations
        for eq in self.ode_equations:
            yd0.append(eq.rhs.eval_with_values(p_vals=par))

        # now append a list of zeros to account for conserved variables
        for eq in self.conservations_variable_list:
            yd0.append(0)

        return yd0

    def construct_fun_DAE(self, t, y0, yd0):

        y_out = np.array([])

        # update dictionary of parameter values to consider initial conditions for Xd_t pools
        for indx, element in enumerate(self.ode_equations):
            self.parameters[str(element.lhs).replace('.', '')] = y0[indx]
            self.parameters[str(element.lhs)] = yd0[indx]

        # update dictionary of parameter values to consider intial conditions for Xd_a_c pools
        for indx2, element in enumerate(self.conservations_variable_list):
            self.parameters[element] = y0[indx + indx2 + 1]

        # start building the y_out vector
        for element in self.ode_equations:
            y_out = np.append(y_out, element.rhs.eval_with_values(p_vals=self.parameters))

        for variable in self.conservations_variable_list:
            eq = self.conservations_equal_zero_dictionary[variable]
            y_out = np.append(y_out,
                              eq.rhs.eval_with_values(p_vals=self.parameters))

        y_out = y_out - yd0

        return y_out

