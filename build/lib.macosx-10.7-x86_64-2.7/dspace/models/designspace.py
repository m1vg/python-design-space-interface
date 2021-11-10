''' Definition of the abstract model class.


'''
import itertools
import sys

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model
from dspace.models.gma import GMASystem
from dspace.models.ssystem import SSystem
from dspace.models.case import Case, CaseIntersection, CaseColocalization
from dspace.models.cyclicalcase import CyclicalCase
from dspace.expressions import Expression
import time
import numpy as np
from functools import cmp_to_key


def sort_cases(x, y):
    x = x.replace('.' , '_')
    y = y.replace('.' , '_')
    x = x.split('_')
    y = y.split('_')
    for i in range(min(len(x), len(y))):
        xi = int(x[i])
        yi = int(y[i])
        if xi < yi:
            return -1
        if xi > yi:
            return 1
    return 0


class DesignSpace(GMASystem):
    
    def __init__(self, equations,
                 parameter_dict=None, 
                 resolve_cycles=False,
                 resolve_instability=False,
                 resolve_conservations=False,
                 number_conservations=0,
                 constraints=None,
                 Xi=None,
                 match_Xi=None,
                 latex_symbols=None,
                 resolve_codominance=False,
                 adjust_codominant_stoichiometry=False,
                 skip_overlapping_codominant_phenotypes=True,
                 consider_mass_balances=False,
                 codominance_adjust_boundaries=False,
                 **kwargs):
        ''' Initializes a new object with the input parameters for a routine
            analysis.
        
        
                 
        Args:
            equations (list): A list of equations in string format defining the 
                system to analyze.
        
        Kwargs:
            parameter_dict (dict): A dictionary of substrings and replacement
                strings to modify original equations.
            
            resolve_cycles (bool): A flag indicating if cycles should be
               automatically resolved. Setting this to true adds significant
               overhead._valid_cases_expand_cyclesDSDesignSpaceCalculateAllValidCases
        '''

        if parameter_dict is not None:
            equations = equations.replace_symbols(parameter_dict)
        super(DesignSpace, self).__init__(equations,
                                          Xi=Xi,
                                          match_Xi=match_Xi,
                                          latex_symbols=latex_symbols, **kwargs)

        setattr(self, '_resolve_cycles', False)
        if constraints is not None:
            if isinstance(constraints, list) is False:
                constraints = [constraints]
            DSDesignSpaceAddConstraints(self._swigwrapper, constraints, len(constraints))
        if resolve_instability is True:
            DSDesignSpaceSetUnstable(self._swigwrapper, True)
        if resolve_conservations is True:
            DSDesignSpaceSetResolveConservations(self._swigwrapper, True)
            DSDesignSpaceSetNumberOfConservations(self._swigwrapper, number_conservations)
        if resolve_codominance is True:
            DSDesignSpaceSetResolveCoDominance(self._swigwrapper, True)
            if adjust_codominant_stoichiometry is True:
                DSDesignSpaceSetAdjustCodominantStoichiometry(self._swigwrapper, True)
            if skip_overlapping_codominant_phenotypes is True:
                DSDesignSpaceSetSkipOverlappingCodominantPhenotypes(self._swigwrapper, True)
        if consider_mass_balances is True:
            fin, fout, signature, metabolic_pools, S_string, rows, columns, rxns = self.generate_mass_balances()
            # print("fin and fout are: ", fin, fout, signature)
            DSDesignSpaceInitializeMassBalances(self._swigwrapper, fin, fout, signature, len(fin),
                                                metabolic_pools._swigwrapper,
                                                S_string, rows, columns, rxns)
        if codominance_adjust_boundaries is True:
            DSDesignSpaceSetAdjustCodominantBoundaries(self._swigwrapper, True)
        if resolve_cycles == True:
            setattr(self, '_resolve_cycles', True)
            DSDesignSpaceCalculateCyclicalCases(self._swigwrapper)
        
    def __del__(self):
        if self._swigwrapper is not None:
            DSDesignSpaceFree(self._swigwrapper)

    def __len__(self):
        return DSDesignSpaceNumberOfCases(self._swigwrapper)+1

    def __getstate__(self):
        odict = self.__dict__.copy()
        odict['_swigwrapper'] = DSSWIGDSDesignSpaceEncodedBytes(self._swigwrapper)
        del odict['_independent_variables']
        return odict

    def __setstate__(self, state):
        self.__dict__.update(state)
        encoded = state['_swigwrapper']
        self.set_swigwrapper(DSSWIGDSDesignSpaceDecodeFromByteArray(encoded))
    
    def _case_with_signature(self, signature, constraints):
        siglist = []
        wild_cards = []
        try:
            case, subcase = signature.split('_')
            signature = case
            subcase = '_' + subcase
        except:
            subcase = ''
            pass
        i = 0
        while i < len(signature):
            if signature[i] == '(':
                start = i+1
                while signature[i] != ')':
                    i += 1
                siglist.append(int(signature[start:i]))
            elif signature[i] == '*':
                num_wild = int(self.signature[i])
                for j in range(num_wild):
                    new_sig = signature.replace('*', str(j+1), 1)
                    wild_cards += self._case_with_signature(new_sig, constraints)
                return wild_cards
            else:
                siglist.append(int(signature[i]))
            i+=1
        index = DSCaseNumberForSignature(siglist, DSDesignSpaceGMASystem(self._swigwrapper))
        return [self(str(index)+subcase, constraints=constraints)]

    def case_number_for_signature(self, signature):
        siglist = []
        i = 0
        if signature == 'Negativevaluefound!':
            return None

        while i < len(signature):
            if signature[i] == '(':
                start = i + 1
                while signature[i] != ')':
                    i += 1
                siglist.append(int(signature[start:i]))
            else:
                siglist.append(int(signature[i]))
            i += 1
        index = DSCaseNumberForSignature(siglist, DSDesignSpaceGMASystem(self._swigwrapper))
        return index
    
    def __call__(self, index_or_iterable, by_signature=False, constraints=None):
        if isinstance(index_or_iterable, (int, str)) is True:
            iterable = [index_or_iterable]
        else:
            iterable = index_or_iterable
        if by_signature is True:
            iterable = [':' + str(i) for i in iterable]
        if constraints is not None:
            if isinstance(constraints, list) is False:
                constraints = [constraints]
        cases = list()
        for index in iterable:
            if isinstance(index, int):
                index = str(index)
            if isinstance(index, str) is True:
                if index[0] == ':':
                    cases += self._case_with_signature(index[1:], constraints)
                    continue
                case_swig = DSDesignSpaceCaseWithCaseIdentifier(self._swigwrapper, index)
                if case_swig is None:
                    raise ValueError('Case "' + index + '" does not exist')
                name = self.name + ': Case ' + index
                case = Case(self,
                            case_swig,
                            name,
                            constraints=constraints)
                eq=Equations(case.equations.system,
                              case.auxiliary_variables)
                cyclical_swig = DSDesignSpaceCyclicalCaseWithCaseIdentifier(self._swigwrapper, index)
                if cyclical_swig is not None:
                    cyclical = CyclicalCase(eq, cyclical_swig,
                                            name=case.name + ' (cyclical)',
                                            latex_symbols=self._latex)
                    cases.append(cyclical)
                else:
                    cases.append(case)
            else:
                raise TypeError('input argument must be a case identifier or case signature')
        if len(cases) == 1:
            if isinstance(index_or_iterable, str) or isinstance(index_or_iterable, int):
                cases = cases[0]
        return cases
            
    def _parse_equations(self, 
                         Xi=None,
                         match_Xi=None, 
                         **kwargs):
        auxiliary_variables = self.auxiliary_variables
        if Xi is not None:
            Xi = [i for i in Xi]
        elif match_Xi is not None:
            Xi = match_Xi.independent_variables
        if Xi is None:
            swigwrapper = DSSWIGDesignSpaceParseWrapper(self.equations.system,
                                                        len(self.equations),
                                                        auxiliary_variables,
                                                        len(auxiliary_variables)
                                                        )
        else:
            swigwrapper = DSSWIGDesignSpaceParseWrapperWithXi(self.equations.system,
                                                              len(self.equations),
                                                              auxiliary_variables,
                                                              len(auxiliary_variables),
                                                              Xi,
                                                              len(Xi),
                                                              )
        self.set_swigwrapper(swigwrapper)
        gma = DSDesignSpaceGMASystem(self._swigwrapper)
        eqs = DSGMASystemEquations(gma)
        equation_list = list()
        for i in range(0, DSGMASystemNumberOfEquations(gma)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equation_list.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        Xda = VariablePool()
        Xda.set_swigwrapper(DSVariablePoolCopy(DSGMASystemXd_a(gma)))
        equations = Equations(equation_list, 
                              auxiliary_variables=Xda.keys(), 
                              latex_symbols=self._latex)
        self._equations = equations

    def update_latex_symbols(self, symbols):
        self._latex.update(symbols)
        gma = DSDesignSpaceGMASystem(self._swigwrapper)
        eqs = DSGMASystemEquations(gma)
        equation_list = list()
        for i in range(0, DSGMASystemNumberOfEquations(gma)):
            expr = DSExpressionAtIndexOfExpressionArray(eqs, i)
            equation_list.append(DSExpressionAsString(expr))
            DSExpressionFree(expr)
        DSSecureFree(eqs)
        Xda = VariablePool()
        Xda.set_swigwrapper(DSVariablePoolCopy(DSGMASystemXd_a(gma)))
        equations = Equations(equation_list, 
                              auxiliary_variables=Xda.keys(), 
                              latex_symbols=self._latex)
        self._equations = equations
                
    def set_swigwrapper(self, ds_swigwrapper):
        self._swigwrapper = ds_swigwrapper
        
        Xd = VariablePool()
        Xd.set_swigwrapper(DSGMASystemXd(DSDesignSpaceGMASystem(ds_swigwrapper)))
        ## for i in VariablePool():
        ##     if i not in self.dependent_variables:
        ##         raise NameError, 'Dependent Variables are inconsistent'
        self._dependent_variables = Xd.copy()
        Xd.set_swigwrapper(None)
        Xi = VariablePool()
        Xi.set_swigwrapper(DSDesignSpaceXi(ds_swigwrapper))
        self._independent_variables = Xi.copy()
        Xi.set_swigwrapper(None)

        Xd_t = VariablePool()
        Xd_t.set_swigwrapper(DSGMASystemXd_t(DSDesignSpaceGMASystem(ds_swigwrapper)))
        self._dependent_variables_no_algebraic = Xd_t.copy()
        Xd_t.set_swigwrapper(None)


    @property
    def dependent_variables(self):
        return self._dependent_variables.keys()

    @property
    def instability(self):
        return DSDesignSpaceUnstable(self._swigwrapper)

    @property
    def dependent_variables_no_algebraic_constraints(self):
        return self._dependent_variables_no_algebraic.keys()

    @property
    def number_of_cases(self):
        return DSDesignSpaceNumberOfCases(self._swigwrapper)

    @property
    def signature(self):
        return DSDesignSpaceSignatureToString(self._swigwrapper)

    @property
    def is_conserved(self):
        return DSDesignSpaceConserved(self._swigwrapper)
    @property
    def number_conservations(self):
        return DSDesignSpaceNumberOfConservations(self._swigwrapper)

    @property
    def _signature(self):
        signature_internal = DSDesignSpaceSignature(self._swigwrapper)
        signature = list()
        for i in range(len(self.equations)*2):
            signature.append(DSUIntegerAtIndexOfIntegerArray(signature_internal, i))
        return signature

    @property
    def should_consider_mass_balances(self):
        return DSDesignSpaceShouldConsiderMassBalances(self._swigwrapper)

    @property
    def number_of_metabolic_blocks(self):
        return DSDesignSpaceNumberOfMetabolicBlocks(self._swigwrapper)

    def fin_at_index(self, i):
        return DSDesignSpaceFinAtIndex(self._swigwrapper, i)

    def fout_at_index(self, i):
        return DSDesignSpaceFoutAtIndex(self._swigwrapper, i)

    def _valid_cases_bounded(self, p_bounds, strict):
        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        for i in lower:
            lower[i] = 1E-20
            upper[i] = 1E20
        for (key,value) in p_bounds.items():
            try:
                min_value,max_value = value
            except TypeError:
                min_value = value
                max_value = value
            if min_value > max_value:
                raise ValueError('parameter slice bounds are inverted: min is larger than max')
            lower[key] = min_value
            upper[key] = max_value
        if strict is True:
            valid_cases = DSDesignSpaceCalculateAllValidCasesForSlice(self._swigwrapper,
                                                                      lower._swigwrapper,
                                                                      upper._swigwrapper)
        else:
            valid_cases = DSDesignSpaceCalculateAllValidCasesForSliceNonStrict(self._swigwrapper,
                                                                               lower._swigwrapper,
                                                                               upper._swigwrapper)
        number_of_cases = DSDictionaryCount(valid_cases)
        cases = list()
        keys = [DSDictionaryKeyAtIndex(valid_cases, i) for i in range(0, number_of_cases)]
        for key in keys:
            case_swigwrapper = DSSWIGVoidAsCase(DSDictionaryValueForName(valid_cases, key))
            cases.append(DSCaseIdentifier(case_swigwrapper))
            DSCaseFree(case_swigwrapper)
        DSDictionaryFree(valid_cases)
        cases.sort(key=cmp_to_key(sort_cases))
        return cases

    def _valid_cases_expanded_bounded(self, p_bounds, strict=True):
        lower = VariablePool(names=self.independent_variables)
        upper = VariablePool(names=self.independent_variables)
        for i in lower:
            lower[i] = 1E-20
            upper[i] = 1E20
        for (key,value) in p_bounds.items():
            try:
                min_value,max_value = value
            except TypeError:
                min_value = value
                max_value = value
            if min_value > max_value:
                raise ValueError('parameter slice bounds are inverted: min is larger than max')
            lower[key] = min_value
            upper[key] = max_value
        valid_cases = DSDesignSpaceCalculateAllValidCasesForSliceByResolvingCyclicalCases(
                       self._swigwrapper,
                       lower._swigwrapper,
                       upper._swigwrapper,
                       strict)
        number_of_cases = DSDictionaryCount(valid_cases)
        cases = list()
        keys = [DSDictionaryKeyAtIndex(valid_cases, i) for i in range(0, number_of_cases)]

        for key in keys:
            case_swigwrapper = DSSWIGVoidAsCase(DSDictionaryValueForName(valid_cases, key))
            cases.append(DSCaseIdentifier(case_swigwrapper))
            DSCaseFree(case_swigwrapper)

        DSDictionaryFree(valid_cases)
        cases.sort(key=cmp_to_key(sort_cases))
        return cases
      
    def _valid_cases_expand_cycles(self, p_bounds, strict=True):
        if p_bounds is not None:
            return self._valid_cases_expanded_bounded(p_bounds, strict=strict)
        valid_cases = DSDesignSpaceCalculateAllValidCasesByResolvingCyclicalCases(self._swigwrapper)
        number_of_cases = DSDictionaryCount(valid_cases)
        cases = list()
        keys = [DSDictionaryKeyAtIndex(valid_cases, i) for i in range(0, number_of_cases)]
        for key in keys:
            case_swigwrapper = DSSWIGVoidAsCase(DSDictionaryValueForName(valid_cases, key))
            cases.append(DSCaseIdentifier(case_swigwrapper))
            DSCaseFree(case_swigwrapper)
        DSDictionaryFree(valid_cases)
        cases.sort(key=cmp_to_key(sort_cases))
        return cases

    def _valid_cases_unstable_subcases(self, cases, case_swigwrapper):

        all_subcases = DSUnstableCaseListAllSubcases(case_swigwrapper, self._swigwrapper)
        number_subcases_valid = DSUnstableCaseSubcasesCount(case_swigwrapper, self._swigwrapper)
        for w in range(0, number_subcases_valid):
            subcase_swigwrapper = DSSWIGVoidAsCase(DSCaseAtIndexOfArray(all_subcases, w))
            #subcase_swigwrapper = DSCaseAtIndexOfArray(all_subcases, w)
            cases.append(DSCaseIdentifier(subcase_swigwrapper))
            DSCaseFree(subcase_swigwrapper)
        if all_subcases is not None:
            DSSecureFree(all_subcases)
        return cases

    def valid_cases(self, p_bounds=None, expand_cycles=True, strict=True):
        if self._resolve_cycles is False:
            expand_cycles = False
        if expand_cycles is True:
            return self._valid_cases_expand_cycles(p_bounds, strict=strict)
        if p_bounds is not None:
            return self._valid_cases_bounded(p_bounds, strict)
        all_cases = DSDesignSpaceCalculateAllValidCases(self._swigwrapper)
        number_valid = DSDesignSpaceNumberOfValidCases(self._swigwrapper) + DSDesignSpaceNumberOfValidBlowingCases(self._swigwrapper, True)
        cases = list()
        for i in range(0, number_valid):
            case_swigwrapper = DSCaseAtIndexOfArray(all_cases, i)
            cases.append(DSCaseIdentifier(case_swigwrapper))
            DSCaseFree(case_swigwrapper)

        if all_cases is not None:
            DSSecureFree(all_cases)
        cases.sort(key=cmp_to_key(sort_cases))
        return cases
    
    def _cyclical_case_as_subcases(self, case_num, case_numbers):
        if case_num not in case_numbers:
            return case_numbers
        case = self(case_num)
        if case.is_cyclical is True:
            case_numbers.remove(case_num)
            new_case_numbers = [str(case_num) + '_' + str(j) for j in range(1, case.number_of_subcases+1)]
            case_numbers = case_numbers + new_case_numbers
            for i in new_case_numbers:
                case_numbers = self._cyclical_case_as_subcases(i, case_numbers)
        return case_numbers
    
    def cycles_to_subcases(self, case_numbers):
        original_cases = case_numbers
        case_numbers = [i for i in case_numbers]
        for i in original_cases:
            case_numbers = self._cyclical_case_as_subcases(i, case_numbers)
        return case_numbers
        
    def valid_intersecting_cases(self, intersects, case_numbers, p_bounds=None, strict=True):
        if isinstance(intersects, list) is False:
            intersects = [intersects]
        if len(case_numbers) == 0:
            return None
        intersections = list()        
        Cases = self(case_numbers)
        ## for i in xrange(len(Cases)):
        ##     case = Cases[i]
        ##     case_num = case.case_number
        ##     case_numbers = self._cyclical_case_as_subcases(case_num, case_numbers)
        valid_cases=set(range(len(case_numbers)))
        if p_bounds is not None:
            lower = VariablePool(names=self.independent_variables)
            upper = VariablePool(names=self.independent_variables)
            for key in lower:
                lower[key] = 1e-20
                upper[key] = 1e20
                for (key,value) in p_bounds.items():
                    try:
                        min_value,max_value = value
                    except TypeError:
                        min_value = value
                        max_value = value
                    if min_value > max_value:
                        raise ValueError('parameter slice bounds are inverted: min is larger than max')
                        lower[key] = min_value
                        upper[key] = max_value
        if 1 in intersects:
            [intersections.append(i) for i in case_numbers if self(i).is_valid(p_bounds=p_bounds, strict=strict) is True]
        sets = [set([i]) for i in valid_cases]
        for i in range(2, max(intersects)+1):
            sets_to_check = sets
            sets = []
            identifiers = range(0, len(sets_to_check))
            comb = itertools.combinations(identifiers, 2)
            for j in comb:
                current_set = set()
                for k in j:
                    current_set = current_set.union(sets_to_check[k])
                if len(current_set) != i:
                    continue
                if current_set in sets:
                    continue
                case_int = CaseIntersection([self(case_numbers[k]) for k in current_set])                
                if case_int.is_valid(p_bounds=p_bounds) is True:
                    if i in intersects:
                        intersections.append([case_numbers[k] for k in current_set])
                    sets.append(current_set)
        return intersections
    
    def co_localize_cases(self, case_numbers, slice_parameters, 
                          constraints=None, p_bounds=None, 
                          optimize=None, minimize=True,
                          project=True,
                          by_signature=False):
        cases = self(case_numbers, by_signature=by_signature)
        to_colocalize = CaseColocalization(cases, slice_parameters, constraints=constraints)
        co_localized = to_colocalize.valid_parameter_set(p_bounds=p_bounds, 
                                                         optimize=optimize,
                                                         minimize=minimize,
                                                         project=project)
        return co_localized
        
        
    def maximum_co_localized_cases(self, case_numbers, slice_variables, p_bounds=None):
        new = False        
        if len(case_numbers) == 0:
            return None
        intersections = list()        
        Cases = self(case_numbers)
        for i in range(len(Cases)):
            case = Cases[i]
            case_num = case.case_number
            case_numbers = self._cyclical_case_as_subcases(case_num, case_numbers)
        valid_cases=range(len(case_numbers))
        if p_bounds is not None:
            lower = VariablePool(names=self.independent_variables)
            upper = VariablePool(names=self.independent_variables)
            for key in lower:
                lower[key] = 1e-20
                upper[key] = 1e20
                for (key,value) in p_bounds.items():
                    try:
                        min_value,max_value = value
                    except TypeError:
                        min_value = value
                        max_value = value
                    if min_value > max_value:
                        raise ValueError('parameter slice bounds are inverted: min is larger than max')
                        lower[key] = min_value
                        upper[key] = max_value
        ## [intersections.append(i) for i in case_numbers if self(i).is_valid()is True]
        sets = [set([i]) for i in valid_cases]
        for i in range(2, len(case_numbers)+1):
            new = True
            sets_to_check = sets
            sets = []
            identifiers = range(0, len(sets_to_check))
            comb = itertools.combinations(identifiers, 2)
            for j in comb:
                current_set = set()
                for k in j:
                    current_set = current_set.union(sets_to_check[k])
                if len(current_set) != i:
                    continue
                if current_set in sets:
                    continue
                case_int = CaseColocalization([self(case_numbers[k]) for k in current_set], slice_variables) 
                pvals=case_int.valid_parameter_set()
                if pvals > 0:
                    if new is True:
                        intersections = []
                    new = False
                    intersections.append([case_numbers[k] for k in current_set])
                    sets.append(current_set)
        return intersections
    
    def intersecting_cases(self, intersects, case_numbers, p_bounds=None, strict=True):
         valid_ints = self.valid_intersecting_cases(intersects, case_numbers, p_bounds=p_bounds, strict=strict)
         if valid_ints is None:
             return None
         case_ints = [CaseIntersection(self(i)) for i in valid_ints]
         return case_ints
    
    def _cyclical_case(self, case, name):
        
        if isinstance(case, int) is False:
            raise TypeError('case must be indicated by its case number')
        sub=DSDesignSpaceCyclicalCaseWithCaseNumber(self._swigwrapper, case)
        if sub is None:
            return None
        case = Case(self, DSDesignSpaceCaseWithCaseNumber(self._swigwrapper, case), name)
        eq6=Equations(case.equations.system, case.auxiliary_variables)
        return CyclicalCase(eq6, sub, name = case.name, latex_symbols=self._latex)
    
    def line_1D_positive_roots(self, function, p_vals, slice_variable, 
                               range_slice, resolution=100):
        
        p_bounds = dict(p_vals)
        p_bounds[slice_variable] = range_slice
        valid_cases = self.valid_cases(p_bounds=p_bounds)
        lines = list()
        X_dict, Y_dict, R_dict = ({},{},{})
        unique_R = set()
        line_styles = ['-', '--', '..']
        colors = ['k', 'r', 'y']
        for case in valid_cases:
            X, Y, R = self(case).line_1D_positive_roots(function, p_vals, slice_variable,
                                                        range_slice, resolution=resolution)
            X_dict[case] = X
            Y_dict[case] = Y
            R_dict[case] = R
            unique_R.update(R)
        for i in X_dict:
            X = X_dict[i]
            Y = Y_dict[i]
            R = R_dict[i]
            for j in unique_R:
                x = list()
                y = list()
                for k in range(len(X)):
                    if R[k] != j:
                        if len(x):
                            lines.append((x, y, j))
                        y=list()
                        x=list()
                        continue
                    x.append(X[k])
                    y.append(Y[k])
                if len(x):
                    lines.append((x, y, j))
        return lines
                
    def data_2D_log_gain_repertoire(self, xaxis, yaxis, zaxis, p_bounds=None, cases=False):
        C=self.valid_cases(p_bounds=p_bounds)
        case_dict = {}
        behavior_set = set()
        stable = list()
        unstable = list()
        for i in C:
            case = self(i)
            key = ''
            p = case.valid_parameter_set()
            x = case.ssystem.log_gain(zaxis, xaxis)
            if x < 0:
                key = '<,'
            elif x == 0:
                key = '0,'
            else:
                key = '>,'
            y = case.ssystem.log_gain(zaxis, yaxis)
            if y < 0:
                key+= '<'
            elif y == 0:
                key += '0'
            else:
                key += '>'
            eigen = case.positive_roots(p)
            if eigen == 0:
                key += ':-'
            else:
                key += ':+'
            if key in case_dict:
                case_dict[key].append(i)
            else:
                case_dict[key] = [i]
            if (x, y, eigen) in behavior_set:
                continue
            behavior_set.add((x, y, eigen))
        if cases is True:
            return case_dict
        return behavior_set

    def dominant_signature(self, xi, results_dic, n):

        # 1. Get keys of results_dic
        keys_dic = results_dic.keys()
        if len(keys_dic) == 0:
            return None

        xd = VariablePool(names=self.dependent_variables)
        dom_sig = []

        # 1.1 If the system has conservations, complement the results_dic with vectors for Xc1, ...Xcn
        if self.number_conservations != 0:
            for c in range(self.number_conservations):
                results_dic.update({'Xc'+str(c+1): np.ones(n)})
                keys_dic = results_dic.keys()

        # 2. Get first element of key and loop over its contents to update the dictionary of dependent variables.
        for ndx in range(n):
            n_dic = {element: results_dic[element][ndx] for element in keys_dic}
            xd.update(n_dic)
            dom_sig.append(DSDesignSpaceDominantSignature(self._swigwrapper, xi._swigwrapper, xd._swigwrapper))

        return dom_sig

    def stoichiometric_connectivity_matrices(self):

        eqs = self._equations
        dependent_variables = self._equations.dependent_variables
        auxiliary_variables = self._equations.auxiliary_variables

        # 1. Get a list of dependent variables excluding auxiliary variables
        metabolic_pools = list(el for el in dependent_variables if el not in auxiliary_variables)

        # 2. Generate a list of unique reaction terms
        rxns = []
        for n, i in enumerate(eqs):
            if str(eqs[n].lhs).replace('.', '') not in metabolic_pools:
                continue
            rhs = Expression(str(eqs[n].rhs))
            for ii in range(DSExpressionNumberOfTerms(rhs._swigwrapper)):
                term = DSExpressionBranchAtIndexAsString(rhs._swigwrapper, ii)
                if term[0] == '-':
                    term = term.replace('-', '', 1)
                if term not in rxns and term != '0':
                    rxns.append(term)

        # 3. initialice the stoichiometric and connectivity matrices
        S = np.zeros([len(metabolic_pools), len(rxns)])
        N = np.zeros([len(metabolic_pools), len(rxns)])

        # 4. Populate the stoichimetric and connectivity matrices
        for n, i in enumerate(eqs):
            pool = str(eqs[n].lhs).replace('.', '')
            if pool not in metabolic_pools:
                continue
            rhs = Expression(str(eqs[n].rhs))
            for ii in range(DSExpressionNumberOfTerms(rhs._swigwrapper)):
                flux = DSExpressionBranchAtIndexAsString(rhs._swigwrapper, ii)
                if flux == '0' or flux == '':
                    continue
                coeff = flux.split('*')[0]
                try:
                    coeff = int(coeff)
                except:
                    coeff = 1 if '-' not in coeff else -1
                if flux[0] == '-':
                    flux = flux.replace('-', '', 1)
                S[metabolic_pools.index(pool), rxns.index(flux)] = + coeff
                N[metabolic_pools.index(pool), rxns.index(flux)] = 1

        return [S, N, metabolic_pools, rxns]

    def process_stoichiometric_matrix(self, S, N, metabolic_pools):

        metabolic_pools = np.array(metabolic_pools)
        for col in range(np.size(S, 1)):
            for row in range(np.size(S, 0)):
                if col == 0:
                    scan = True if S[row, col] == 0.0 else False
                else:
                    scan = True if S[row, col] == 0.0 and sum(np.abs(S[row, 0:col])) == 0.0 else False
                if scan is True:
                    for row_swap in range(row, np.size(S, 0)):
                        if S[row_swap, col] != 0.0:
                            S[[row, row_swap]] = S[[row_swap, row]]
                            N[[row, row_swap]] = N[[row_swap, row]]
                            metabolic_pools[[row, row_swap]] = metabolic_pools[[row_swap, row]]

        return S, N, list(metabolic_pools)

    def generate_mass_balances(self):

        # 1. We first get the stoichiometric and connectivity matrices and then process the stoichiometric matrix
        S, N, metabolic_pools, rxns = self.stoichiometric_connectivity_matrices()
        S, N, metabolic_pools = self.process_stoichiometric_matrix(S, N, metabolic_pools)

        # print("The stoichiometric matrix is: ", S)
        # print("metabolic_pools: ", metabolic_pools)
        # print("rxns: ", rxns)

        # 2. new: we now identify metabolic groups
        assigned = []
        groups = []
        for n, pool_i in enumerate(metabolic_pools):
            index1 = metabolic_pools.index(pool_i)
            # the first element is used to create the first group
            if n == 0:
                group_i = [index1]
                groups.append(group_i)
                assigned = [index1]
                continue
            if pool_i in assigned:
                continue

            # Now we check if pool_i is connected to any group_i in groups
            new_group = True
            for group_i in groups:
                # print("checking if pool {} is connected with the groups: ".format(pool_i), groups)
                for index_assigned in group_i:
                    # print("checking index_assigned: ", index_assigned)
                    # print("building the product: ",  N[index1, :], N[index_assigned, :], "=", sum(N[index1, :] * N[index_assigned, :]))
                    if sum(N[index1, :] * N[index_assigned, :]) != 0.0:
                        group_i.append(index1)
                        assigned.append(index1)
                        new_group = False
                        # print("should not create a new group for this pool")
                        break
                    else:
                        continue
            if new_group is True:
                group_i = [index1]
                groups.append(group_i)
                assigned.append(index1)

        # 3. we now loop for each group and indentify the mass balances
        fin_vector = []
        fout_vector = []
        metabolic_pools_dic = {}
        for pool_n, group_i in enumerate(groups):
            S_i = S[group_i, :]
            N_i = N[group_i, :]
            metabolic_pools_i = [metabolic_pools[i] for i in group_i]
            metabolic_pools_dic.update(dict((met, pool_n) for met in metabolic_pools_i))
            fin, fout = self.indentify_input_output_fluxes(S_i, N_i, rxns, metabolic_pools_i)
            fin_vector.append(fin)
            fout_vector.append(fout)
        metabolic_pools_variable_pool = VariablePool()
        metabolic_pools_variable_pool.update(metabolic_pools_dic)

        # If auxiliary variables are present, include those with the pool of the first rxn in which they appear.
        self.define_pools_for_auxiliary_variables(metabolic_pools_variable_pool, rxns, groups, S)

        fin_string = []
        fout_string = []
        signature = []
        for n, fin_i in enumerate(fin_vector):
            fin_string.append('Fin = ' + '+'.join(fin_i))
            signature.append(str(len(fin_i)))
            signature.append('0')
        for n, fout_i in enumerate(fout_vector):
            fout_string.append('Fout = ' + '+'.join(fout_i))
            signature[2*n+1] = str(len(fout_i))

        # The stoichiometric matrix is recalculated to match the order of the metabolites
        S, N, metabolic_pools, rxns = self.stoichiometric_connectivity_matrices()
        S_string = []
        rows = np.size(S, 0)
        columns = np.size(S, 1)

        for row in range(rows):
            for col in range(columns):
                S_string.append(str(S[row, col]))

        return fin_string, fout_string, signature, metabolic_pools_variable_pool, S_string, rows, columns, rxns

    def indentify_input_output_fluxes(self, S, N, rxns, metabolic_pools):

        # Exchange reactions are identified
        sum_rows = np.sum(N, axis=0)
        exchange_rxns = np.array(rxns)[sum_rows == 1.0]

        # we now loop exchange_rxns and construct fin and fout vectors
        f_in = []
        f_out = []
        for rxn in exchange_rxns:
            column = rxns.index(rxn)
            for n in range(len(metabolic_pools)):
                if (np.sign(S[n, column]) == -1.0):
                    f_out.append(rxn)
                    break
                if (np.sign(S[n, column]) == 1.0):
                    f_in.append(rxn)
                    break

        # we now loop over the stoichiometric matrix looking for exchange reaction involving multiple substrates.
        for column, rxn in enumerate(rxns):
            pos = 0
            neg = 0
            for row in range(len(metabolic_pools)):
                if S[row, column] == 0.0:
                    continue
                if np.sign(S[row, column]) == 1.0:
                    pos += 1
                else:
                    neg += 1
            if pos == 0 and neg != 0 and rxn not in f_out:
                f_out.append(rxn)
            if pos != 0 and neg == 0 and rxn not in f_in:
                f_in.append(rxn)

        return f_in, f_out

    def define_pools_for_auxiliary_variables(self, metabolic_pools_variable_pool, rxns, groups, S):

        # we need to assign metabolic_pools to the reactions.
        # print("The metabolic pools are: ", metabolic_pools_variable_pool)
        metabolic_pools_rxns = {}
        rxns = np.array(rxns)
        for n, group_i in enumerate(groups):
            S_i = S[group_i, :]
            mask = np.sum(abs(S_i), 0) != 0
            rxns_block_i = rxns[mask]
            for rxns_i in rxns_block_i:
                if rxns_i in metabolic_pools_rxns.keys():
                    er = "Reaction belongs to multiple metabolic pools: " + str(rxns_i)
                    # print("The metabolic pools are: ", metabolic_pools_variable_pool)
                    raise ValueError(er)
                else:
                    metabolic_pools_rxns.update({rxns_i: n})
        aux = self.auxiliary_variables
        for eq in self.equations:
            # get the left hand side of each equation
            var = str(eq.lhs).replace('.', '')
            n_terms = len(str(eq.rhs).split('+'))
            if var in aux and n_terms > 1:
                # Then we need to complement metabolic_pools_variable_pool with auxiliary variables and their pools
                if var in metabolic_pools_variable_pool.keys():
                    raise ValueError("Auxiliary Variable belongs to multiple metabolic pools")
                else:
                    for rxns_i in metabolic_pools_rxns.keys():
                        if var in rxns_i:
                            metabolic_pools_variable_pool.update({var: metabolic_pools_rxns[rxns_i]})


