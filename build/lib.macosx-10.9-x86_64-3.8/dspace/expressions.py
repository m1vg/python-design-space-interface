from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool

class Expression(object):
    
    def __init__(self, string_repr):
        
        setattr(self, '_swigwrapper', None)
        if string_repr is not None:
            self._swigwrapper = DSExpressionByParsingString(string_repr)
    
    def __del__(self):
        
        if self._swigwrapper is not None:
            DSExpressionFree(self._swigwrapper)
        
    def __str__(self):
        
        return DSExpressionAsString(self._swigwrapper)
        
    def latex(self, substitution_dictionary=None):
        
        subs=dict(log=r'\log',
                  log10=r'\log_{10}',
                  ln=r'\ln')
        if substitution_dictionary is not None:
            subs.update(substitution_dictionary)
        latex_dict = DSSWIGDSDictionaryFromPyDict(subs)
        string = DSExpressionAsLatexString(self._swigwrapper, latex_dict)
        DSSWIGDSDictionaryFreeCharValues(latex_dict)
        return string
        
    def __repr__(self):
        
        return str(self)
        
    def subst(self, p_vals=None, **kwargs):
        
        if p_vals is None:
            p_vals = VariablePool(kwargs)
        expr = Expression(None)
        expr._swigwrapper = DSExpressionByCompressingConstantVariables(self._swigwrapper,
                                                                       p_vals._swigwrapper)
        return expr
    
    def subst(self, replace_dict, **kwargs):
        replace_dict = dict(replace_dict)
        replace_dict.update(kwargs)
        new_expression = Expression('0')
        DSExpressionFree(new_expression._swigwrapper)
        new_expression._swigwrapper = DSExpressionCopy(self._swigwrapper)
        for key,value in replace_dict.items():
            target = Expression(str(key))
            subs = Expression(str(value))
            temp = DSExpressionByReplacingSubExpression(
                    new_expression._swigwrapper,
                    target._swigwrapper,
                    subs._swigwrapper)
            DSExpressionFree(new_expression._swigwrapper)
            new_expression._swigwrapper = temp
        return new_expression
        
    def eval_with_values(self, p_vals=None, **kwargs):
        
        if p_vals is None:
            p_vals = VariablePool(kwargs)
        return DSExpressionEvaluateWithVariablePool(self._swigwrapper,
                                                    p_vals._swigwrapper)
    
    @property
    def lhs(self):
        expr = None
        lhs_swigwrapper = DSExpressionEquationLHSExpression(self._swigwrapper)
        if lhs_swigwrapper is not None:
            expr = Expression(None)
            expr._swigwrapper = lhs_swigwrapper
        return expr
    
    @property  
    def rhs(self):
        expr = None
        rhs_swigwrapper = DSExpressionEquationRHSExpression(self._swigwrapper)
        if rhs_swigwrapper is not None:
            expr = Expression(None)
            expr._swigwrapper = rhs_swigwrapper
        return expr
        
    @property
    def variables(self):
        vp = DSExpressionVariablesInExpression(self._swigwrapper)
        pv = VariablePool()
        pv.set_swigwrapper(vp)
        return pv.keys()

    