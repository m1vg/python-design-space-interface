''' Definition of the abstract model class.


'''
import itertools

from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.models.designspace import DesignSpace
import numpy as np

class GraphGenerator(object):
    
    def __init__(self, design_space, included_variables=[]):
        
        setattr(self, '_design_space', design_space)
        setattr(self, '_included_variables', included_variables)
        setattr(self, '_flux_nodes', dict())
        variables = list()
        for i in design_space.dependent_variables:
            variables.append(i)
        for i in design_space.independent_variables:
            variables.append(i)
        setattr(self, '_variables', variables)
        
    def graph_description(self, included_independent_variables=None):
        ds= self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)
        connectivity = DSGMASystemNetworkConnectivity(gma)
        influences = list()
        gv_string = ''
        signature = ds._signature
        for i in xrange(len(connectivity)):
            k = i
            term = connectivity[i]
            for j in xrange(len(signature)):
                if (k // signature[j]) > 0:
                    k -= signature[j]
                else:
                    k = k % signature[j]
                    key = (j, k)
                    break
            influences.append((key,[j for j in xrange(len(term)) if term[j] > 0])) 
        flux_data = dict()
        for i in xrange(len(influences)):
            if influences[i][0][0] % 2 == 1:
                data = [self._variables[influences[i][0][0]//2], '->']
            else:
                data = ['->', self._variables[influences[i][0][0]//2]]
            influence = list()
            for j in influences[i][1]:
                if j >= len(self._design_space.equations): 
                    if self._variables[j] not in self._included_variables:
                       continue
                influence.append(self._variables[j])
            flux_data[i] = [[data], influence]
        pp_sets = list()
        for i in xrange(len(self._design_space.equations)):
            for j in xrange(i, len(self._design_space.equations)):
                ppr1 = DSGMASystemPrecursorProductRelationships(gma, i, j)
                ppr2 = DSGMASystemPrecursorProductRelationships(gma, j, i)
                if ppr1 is None:
                    ppr1 = []
                if ppr2 is None:
                    ppr2 = []
                ppr1 = zip(*ppr1)
                ppr2 = zip(*ppr2)
                for p1,p2 in ppr1+ppr2:
                    in_set = False
                    p1 = int(p1)
                    p2 = int(p2)
                    for a_set in pp_sets:
                        if len(a_set.intersection([p1,p2])) > 0:
                            a_set.update([p1,p2])
                            in_set = True
                            break
                    if in_set is False:
                        pp_sets.append(set([p1,p2]))
        for i in pp_sets:
            first_index = None
            for j in i:
                if first_index is None:
                    first_index = j
                    continue
                flux_data[first_index][0] += flux_data[j][0]
                flux_data[first_index][1] += flux_data[j][1]
                flux_data.pop(j)
        print 'digraph {\n    graph[layout=neato,normalize=true];\n    node[shape=plaintext];\n    edge[weight=2]'
        for i in flux_data:
            print '    ' + str(i) + '[shape=circle,width=.0,height=.0,label=""];'
            data = flux_data[i][0]
            external = set(flux_data[i][1])
            for j in data:
                ## print external, j[0], j[1]
                if j[0] == '->':
                    if len(data) == 1:
                        print '    start' + str(i) + ' -> ' + str(i) + '[arrowhead=none];'
                        print '    start' + str(i) + ' [shape=circle,width=.01,height=.01,label=""];'
                    print '    ' + str(i) + j[0] + j[1] + ';'
                    if external.issuperset([j[1]]):
                        external.remove(j[1])
                else:
                    print '    ' + j[0] + j[1] + str(i) + '[arrowhead=none];'
                    if len(data) == 1:
                        print '    ' + str(i) + ' -> end' + str(i) + ';'
                        print '     end' + str(i)+ ' [shape=circle,width=.01,height=.01,label=""];'
                    if external.issuperset([j[0]]):
                        external.remove(j[0])
            for k in external:
                if k in self._design_space.dependent_variables:
                    print '    ' + k + '->' + str(i) + '[weight=1];'
                elif k in self._included_variables:
                    print '    ' + k + '->' + str(i) + '[weight=1];'                        
        print '};'
        return flux_data
        
    def connectivity(self):
        ds= self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)
        connectivity = DSGMASystemNetworkConnectivity(gma)
        return connectivity
    
    def all_fluxes(self):
        ds= self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)
        flux_dict = DSGMASystemFluxDictionary(gma)
        number_of_fluxes = DSDictionaryCount(flux_dict)
        fluxes = dict()
        keys = [DSDictionaryKeyAtIndex(flux_dict, i) for i in xrange(0, number_of_fluxes)]
        for key in keys:
            expr = DSSWIGVoidAsExpression(DSDictionaryValueForName(flux_dict, key))
            flux = DSExpressionAsString(expr)
            DSExpressionFree(expr)
            flux = flux.strip('-')
            fluxes[int(key)] = flux
        DSDictionaryFree(flux_dict)
        return fluxes