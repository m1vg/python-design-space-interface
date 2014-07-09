''' Definition of the abstract model class.


'''

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.base import Equations,Model



class GMASystem(Model):
    
    def __init__(self, equations, name=None, swigwrapper=None, **kwargs):
        super(GMASystem, self).__init__(equations, name=name, **kwargs)
        setattr(self, '_swigwrapper', None)
        setattr(self, '_independent_variables', None)

        if swigwrapper is not None:
            self.set_swigwrapper(swigwrapper)
        else:
            self._parse_equations(**kwargs)
        
    def __del__(self):
        if self._swigwrapper is not None:
            DSGMASystemFree(self._swigwrapper)
        
    def set_swigwrapper(self, gma_swigwrapper):
        self._swigwrapper = gma_swigwrapper
        
        Xd = VariablePool()
        Xd.set_swigwrapper(DSGMASystemXd(gma_swigwrapper))
        for i in VariablePool():
            if i not in self.dependent_variables:
                raise NameError, 'Dependent Variables are inconsistent'

        Xi = VariablePool()
        Xi.set_swigwrapper(DSVariablePoolCopy(DSGMASystemXi(gma_swigwrapper)))
        self._dependent_variables = Xd.copy()
        self._independent_variables = Xi
        Xd.set_swigwrapper(None)
        
    def _parse_equations(self, **kwargs):
        auxiliary_variables = self.auxiliary_variables
        swigwrapper = DSSWIGGMASystemParseWrapper(self.equations,
                                                  len(self.equations),
                                                  auxiliary_variables,
                                                  len(auxiliary_variables)
                                                  )
        self.set_swigwrapper(swigwrapper)
    
    @property
    def independent_variables(self):
        return self._independent_variables.keys()
      