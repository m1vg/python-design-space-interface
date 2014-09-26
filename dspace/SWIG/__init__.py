from distutils.version import StrictVersion
from dspace.SWIG.dspace_interface import DSDesignSpaceToolboxVersionString

compatability = ['0.2.0a1', '-']
__version__ = DSDesignSpaceToolboxVersionString()

if StrictVersion(__version__) < StrictVersion(compatability[0]):
    raise ImportError, 'Needs Design Space Toolbox V2 C Library >= '+compatability[0]
    
if len(compatability[1]) > 1:
    if StrictVersion(__version__) > StrictVersion(compatability[1]):
        raise ImportError, 'Needs Design Space Toolbox V2 C Library <= '+compatability[1]
    
