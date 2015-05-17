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
        variables += design_space.dependent_variables + design_space.independent_variables
        setattr(self, '_variables', variables)
        
    def flux_identifiers(self):
        ds = self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)
        connectivity = self.connectivity()
        fluxes = self.all_fluxes()
        flux_identifiers = {i:i for i in fluxes}
        equivalence = DSGMASystemEquivalentFluxes(gma)
        if equivalence is None:
            return flux_identifiers
        for i in xrange(len(equivalence)):
            for j in xrange(len(equivalence[i])):
                if equivalence[i][j] != 0.0:
                    for k in xrange(i):
                        if equivalence[k][j] != 0.0:
                            flux_identifiers[i] = k
                            break
                    break
        return flux_identifiers
        
    def network_data(self):
        ds = self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)       
        connectivity = self.connectivity()
        flux_identifiers = self.flux_identifiers()
        network_data = {}
        variables = {}
        index = 0
        internal_arrow = '[arrowhead=none]'
        for i,variable in enumerate(ds.dependent_variables):
            if variable in ds.auxiliary_variables:
                continue
            for positive in xrange(ds._signature[2*i]):
                key = flux_identifiers[index+positive]
                link = str(key) + ' -> ' + variable
                if key in network_data:
                    network_data[key]['positive'].append(link)
                else:
                    variables[key] = []
                    network_data[key] = {'positive':[link],
                                         'negative':[]}
            index += ds._signature[2*i]
            for negative in xrange(ds._signature[2*i+1]):
                key = flux_identifiers[index+negative]
                link = variable + ' -> ' + str(key)+internal_arrow
                if key in network_data:
                    variables[key].append(variable)
                    network_data[key]['negative'].append(link)
                else:
                    variables[key] = [variable]
                    network_data[key] = {'positive':[],
                                         'negative':[link]}
            index += ds._signature[2*i+1]
        for key in network_data:
            if len(network_data[key]['negative']) == 0:
                network_data[key]['negative'].append('start_'+str(key)+' -> ' + str(key)+internal_arrow)
            if len(network_data[key]['positive']) == 0:
                network_data[key]['positive'].append(str(key)+' -> end_'+ str(key))
            network_data[key] = network_data[key]['positive'] + network_data[key]['negative']
        return network_data, variables
        
    def graph_regulation(self, included_variables):
        ds = self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)       
        all_variables = ds.dependent_variables + ds.independent_variables
        show_variables = [i for i in ds.dependent_variables if i not in ds.auxiliary_variables]
        show_variables += included_variables
        network_data,variable_links = self.network_data()
        connectivity = self.connectivity()
        flux_identifiers = self.flux_identifiers()
        regulation = []
        positive_r = '[arrowhead=vee]'
        negative_r = '[arrowhead=onormal]'
        for flux in xrange(len(flux_identifiers)):
            key = flux_identifiers[flux]
            for j in xrange(len(connectivity[flux])):
                variable = all_variables[j]
                if variable in variable_links[key]:
                    continue
                if variable not in show_variables:
                    continue
                if connectivity[flux][j] > 0.0:
                    link = variable + ' -> ' + str(key) + positive_r
                elif connectivity[flux][j] < 0.0:
                    link = variable + ' -> ' + str(key) + negative_r
                else:
                    continue
                if link not in regulation:
                    regulation.append(link)
        return regulation            
        
    def graph_properties(self, graph_type):
        properties =  'graph[layout='+graph_type + ',normalize=true];'
        properties += 'node[shape=plaintext];'
        ## properties += 'rankdir=LR;'
        return properties
        
    def subgraph_properties(self, key):
        properties = 'rank=same;'
        properties += 'color=none;'
        properties += 'edge[weight=10];'
        terminal_node_p = '[style=invis,shape=point,label=""];'#shape=circle,width=.01,height=.01,label=""];'
        inner_node_p = '[shape=circle,width=.01,height=.01,label=""];'
        properties += 'start_'+str(key) + terminal_node_p
        properties += 'end_'+str(key) + terminal_node_p
        properties += str(key) + inner_node_p
        return properties
        
    def graph_description(self, graph_type='dot', included_variables=[]):
        network_data,variable_links = self.network_data()
        regulation = self.graph_regulation(included_variables)
        graph_string = 'digraph {'
        graph_string += self.graph_properties(graph_type)
        for key in network_data:
            graph_string += 'subgraph cluster_'+ str(key) + ' {'
            graph_string += self.subgraph_properties(key)
            graph_string += ';'.join(network_data[key])
            graph_string += '}'
        graph_string += 'subgraph {rankdir=LR;'#constraint=false;concentrate=true;'
        graph_string += ';'.join(regulation)
        graph_string += '}}'
        return graph_string
               
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