from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.expressions import Expression
import dspace.models.base

from dspace.plotutils.monkey_patching import monkeypatch_method
from IPython.display import display, Math

@monkeypatch_method(dspace.models.base.Equations)
def _repr_latex_(self):
    string = r'\begin{array}{l}'
    for i in self.system:
        eq = Expression(i)
        string += eq.__latex_str__(substitution_dictionary=self._latex)
        string += r'\\'
    string += r'\end{array}'
    return string
