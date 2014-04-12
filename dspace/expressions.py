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
    
    