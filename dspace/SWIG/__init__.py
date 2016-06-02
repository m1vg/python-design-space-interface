from distutils.version import StrictVersion

from dspace.SWIG.dspace_interface import DSSWIGAssignErrorFunctions
from dspace.SWIG.dspace_interface import DSDesignSpaceToolboxVersionString

__supported_versions__ = ['0.3.0a6', '-']
__version__ = DSDesignSpaceToolboxVersionString()

if StrictVersion(__version__) < StrictVersion(__supported_versions__[0]):
    raise ImportError, 'Needs Design Space Toolbox V2 C Library >= '+__supported_versions__[0]
    
if len(__supported_versions__[1]) > 1:
    if StrictVersion(__version__) >= StrictVersion(__supported_versions__[1]):
        raise ImportError, 'Needs Design Space Toolbox V2 C Library < ' +  __supported_versions__[1]    

DSSWIGAssignErrorFunctions()
