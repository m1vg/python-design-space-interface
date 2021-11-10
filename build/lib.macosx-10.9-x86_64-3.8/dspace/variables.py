''' A higher level interface to the DSVariablePool C construct. This
    interface uses the SWIG generated wrapper for the designspace C.
    The variable pool subclasses the dict object, but allows access
    to the variable pool by either key or index, as both are useful
    in the C library. '''

from collections import OrderedDict

##
# Required C api functionality.
##
SWIG_REQUIREMENTS = ['DSVariablePoolSetReadWriteAdd',
                     'DSVariablePoolFree',
                     'DSVariablePoolNumberOfVariables',
                     'DSVariablePoolVariableAtIndex',
                     'DSVariablePoolAlloc',
                     'DSVariablePoolSetValueForVariableWithName',
                     'DSVariablePoolAddVariableWithName',
                     'DSVariablePoolHasVariableWithName',
                     'DSVariablePoolIndexOfVariableWithName'
                     ]

module = __import__('dspace.SWIG.dspace_interface', fromlist=SWIG_REQUIREMENTS)

for function in SWIG_REQUIREMENTS:
    globals()[function] = getattr(module,function)

class VariablePool(dict):
    ''' A python class that serves as a wrapper to the DSVariablePool object
        in the designspace C library. The Variable Pool object is a subclass of
        dict and is used to reference dependent and independent variables, as 
        well as system parameters. The dict object is ordered by the internal
        C data structure.
    '''
    
    def __init__(self, names=None, **kwargs):
        ''' Init method only accepts a list of names. A dictionary cannot be used
            to initialize a VariablePool because dictionaries are unordered. The
            order of the variables in a VariablePool is important and so initialization
            is restricted by the use of iterables of the form of key value pairs.
        '''
        super(VariablePool, self).__init__()
        setattr(self, '_swigwrapper', None)
        setattr(self, '_keys', list())
        if isinstance(names, list) is True:
            names = OrderedDict([(key,1.) for key in names])
        self.update(names=names, **kwargs)
            
    def update(self, names=None, **kwargs):        
        if names is not None:
            for key, value in names.items():
                self[key] = value
        for key, value in kwargs.items():
            self[key] = value
    
    def set_swigwrapper(self, swigwrapper):
        self._swigwrapper = swigwrapper
        if self._swigwrapper is None:
            return
        for i in range(0, DSVariablePoolNumberOfVariables(swigwrapper)):
                variable = DSVariablePoolVariableAtIndex(swigwrapper, i)
                self[variable[0]] = variable[1]

    def __del__(self):
        ''' Method to destroy a VariablePool object. The wrapped C construct must be
            freed to avoid memory leaks.'''
        if self._swigwrapper is None:
            return
        DSVariablePoolSetReadWriteAdd(self._swigwrapper)
        DSVariablePoolFree(self._swigwrapper)
        
    def __setattr__(self, name, value):
        ''' Restricts the attribute modification. The VariablePool can only have a
            _swigwrapper attribute, which is immutable once assigned. '''
        super(VariablePool, self).__setattr__(name, value)
    
    def __getstate__(self):
        odict = [(i,self[i]) for i in self.keys()]
        return odict
    
    def __setstate__(self, state):
        setattr(self, '_swigwrapper', None)
        for key, value in state:
            self[key] = value
        
    def __setitem__(self, name, value):
        if hasattr(self, '_swigwrapper') is False:
            setattr(self, '_swigwrapper', None)
        if hasattr(self, '_keys') is False:
            setattr(self, '_keys', list())
        if self._swigwrapper == None:
            self._swigwrapper = DSVariablePoolAlloc()
        if isinstance(name, str) is False:
            raise TypeError('VariablePool keys must be strings')
        if DSVariablePoolHasVariableWithName(self._swigwrapper, name) == False:
            DSVariablePoolAddVariableWithName(self._swigwrapper, name)
        if name not in self._keys:
            self._keys.append(name)
        value = float(value)
        DSVariablePoolSetValueForVariableWithName(self._swigwrapper,
                                                  name,
                                                  value)
        super(VariablePool, self).__setitem__(name, value)
    
    def keys(self):
        keys = list()
        return keys+[i for i in self._keys]
    
    def iterkeys(self):
        return iter(i for i in self.keys())
        
    def viewkeys(self):
        return self.iterkeys()
                    
    def copy(self):
        newPool = VariablePool()
        for i in self._keys:
            newPool[i] = self[i]
        return newPool

    def index(self, name):
        return self.keys().index(name)
        
