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
        
    def __latex_str__(self, substitution_dictionary=None):
        
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