from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.expressions import Expression
import dspace.models.base

from dspace.plotutils.monkey_patching import monkeypatch_method
from IPython.display import display, Math

@monkeypatch_method(dspace.models.base.Equations)
def _repr_latex_(self):
    if len(self) == 1:
        eq = self._eq[0]
        string = '\[' + eq.latex(substitution_dictionary=self._latex) + '\]'
        return string
    string = r'\begin{equation}'
    for eq in self._eq:
        string += eq.latex(substitution_dictionary=self._latex)
        string += r'\\'
    string += r'\end{equation}'
    return string
