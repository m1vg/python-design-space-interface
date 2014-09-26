''' Definition of abstract base class 'Model' and the Equations object container
    for systems of equations.
    


'''
import string

SWIG_REQUIREMENTS = ['DSSSystemNumberOfEquations',
                     'DSSSystemSolution',
                     'DSSSystemEquations',
                     'DSSSystemXi',
                     'DSSSystemXd',
                     'DSSSystemLogarithmicGain',
                     'DSSWIGSSystemParseWrapper',
                     'DSGMASystemXi',
                     'DSGMASystemXd',
                     'DSSWIGGMASystemParseWrapper',
                     'DSExpressionAtIndexOfExpressionArray',
                     'DSExpressionAsString'
                     ]

module = __import__('dspace.SWIG.dspace_interface', fromlist=SWIG_REQUIREMENTS)

for function in SWIG_REQUIREMENTS:
    globals()[function] = getattr(module,function)


from dspace.variables import VariablePool

class Equations(object):
    
    def __init__(self, system, auxiliary_variables=[], latex_symbols=None):
        ''' Init method for the model base class.
        
        The model object is initialized with data to construct
        a dictionary of key value pairs, where the keys represent
        names of the model's dependent variables, and the values
        are the equations in the system.
        '''
        setattr(self, '_eq', list())
        setattr(self, '_auxiliary_variables', list())
        setattr(self, '_latex', dict())
        if isinstance(system, list) is False:
            system = [system]
        if isinstance(auxiliary_variables, list) is False:
            auxiliary_variables = [auxiliary_variables]            
        for i in system:
            if isinstance(i, str) is False:
                raise TypeError, 'ODE must be a string'
            self._eq.append(i)
        for i in auxiliary_variables:
            if isinstance(i, str) is False:
                raise TypeError, 'ODE must be a string'
            self._auxiliary_variables.append(i)
        if latex_symbols is not None:
            self._latex.update(latex_symbols)
    
    @property                    
    def system(self):
        return list(self._eq)

    @property
    def dependent_variables(self):
        differentiated = list()
        for i in self._eq:
            temp = i.split('=')[0].strip()
            temp = temp.split('.')
            if len(temp) == 1 or temp[0] not in (string.lowercase + string.uppercase + '_') is True:
                continue
            differentiated.append(temp[0])
        differentiated += self._auxiliary_variables
        return differentiated
    
    @property
    def auxiliary_variables(self):
        return list(self._auxiliary_variables)
    
    def __getitem__(self, index):
        return self._eq[index]
    
    def __len__(self):
        return len(self._eq)
        
    def __repr__(self):
        return str(self._eq)
        
    def replace_symbols(self, symbol_dict):
        
        eq = [i for i in self._eq]
        for i in xrange(0, len(self._eq)):
            for j in symbol_dict:
                eq[i] = eq[i].replace(j, str(symbol_dict[j]))
        eqs = Equations(eq, self._auxiliary_variables)
        return eqs
    

class Model(object):
    
    def __init__(self, equations, name=None, description=None, latex_symbols=None, **kwargs):
        ''' Init method for the model base class.
        
        The model object is initialized with data to construct
        a dictionary of key value pairs, where the keys represent
        names of the model's dependent variables, and the values
        are the equations in the system.
        '''
        setattr(self, '_equations', equations)
        if name is None:
            name = 'Unnamed'
        setattr(self, '_name', name)
        if description is None:
            description = ''
        setattr(self, '_description', description)
        setattr(self, '_latex', dict())
        if latex_symbols is not None:
            self._latex.update(latex_symbols)
        else:
            self._latex.update(equations._latex)


    @property
    def name(self):
        return self._name
    
    @property
    def equations(self):
        return self._equations
    
    @property
    def auxiliary_variables(self):
        return self._equations.auxiliary_variables
    
    @property
    def dependent_variables(self):
        return self._equations.dependent_variables
    
    def __repr__(self):
        return 'Model: ' + self.name + '\n' + str(self.equations)     
   
            