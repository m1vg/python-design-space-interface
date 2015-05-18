''' Definition of the abstract model class.


'''
import itertools


from dspace.SWIG.dspace_interface import *
from dspace.variables import VariablePool
from dspace.expressions import Expression
from dspace.models.designspace import DesignSpace
import numpy as np
from math import *

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
    
    def flux_concentrations(self, p_vals):
        ds = self._design_space
        case = self._design_space.valid_cases(p_bounds=p_vals)[0]
        case = self._design_space(case)
        ss = case.steady_state(p_vals)
        values = VariablePool()
        values.update(ss)
        values.update(p_vals)
        fluxes = self.all_fluxes()
        max_value = None
        min_value = None
        flux_value = {}
        for i in fluxes:
            ex = Expression(fluxes[i])
            flux_value[i] = log10(ex.eval_with_values(values))
            if max_value is None:
                max_value = flux_value[i]
            else:
                max_value = max(flux_value[i], max_value)
            if min_value is None:
                min_value = flux_value[i]
            else:
                min_value = min(flux_value[i], min_value)
        return flux_value, [min_value, max_value]
        
    def network_data(self, p_vals, cmap=None):
        ds = self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)       
        connectivity = self.connectivity()
        flux_identifiers = self.flux_identifiers()
        network_data = {}
        variables = {}
        index = 0
        external_arrow = '[color={0}]'
        internal_arrow = '[arrowhead=none,color={0}]'
        terminal_node_p = '[style=invis,shape=point,label=""]'
        inner_node_p = '[shape=circle,width=.01,height=.01,label="",color={0}]'
        if p_vals is not None and cmap is not None:
            concentrations, [vmin, vmax] = self.flux_concentrations(p_vals)
            colors={}
            for i in flux_identifiers:
                key = flux_identifiers[i]
                color = cmap((concentrations[key]-vmin)/(vmax-vmin), bytes=True)
                color_sets = [(hex(j).split('x')[1]*2)[:2] for j in color]
                color = '"#'+''.join(color_sets)+'"'
                colors[key] = color
        else:
            colors = {flux_identifiers[i]:'"black"' for i in flux_identifiers}
        for i,variable in enumerate(ds.dependent_variables):
            if variable in ds.auxiliary_variables:
                continue
            for positive in xrange(ds._signature[2*i]):
                key = flux_identifiers[index+positive]
                link = str(key) + ' -> ' + variable + external_arrow.format(colors[key])
                if key in network_data:
                    network_data[key]['positive'].append(link)
                else:
                    variables[key] = []
                    network_data[key] = {'positive':[link],
                                         'negative':[]}
            index += ds._signature[2*i]
            for negative in xrange(ds._signature[2*i+1]):
                key = flux_identifiers[index+negative]
                link = variable + ' -> ' + str(key)+internal_arrow.format(colors[key])
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
                node = 'start_'+str(key)+terminal_node_p
                link = 'start_'+str(key)+' -> ' + str(key)+internal_arrow.format(colors[key])
                network_data[key]['negative'] += [link, node]
            if len(network_data[key]['positive']) == 0:
                node = 'end_'+str(key)+terminal_node_p
                link = str(key)+' -> end_'+ str(key)+external_arrow.format(colors[key])
                network_data[key]['positive'] += [link, node]
            node = str(key) + inner_node_p.format(colors[key])
            network_data[key]['positive'].append(node)
            network_data[key] = network_data[key]['positive'] + network_data[key]['negative']
        return network_data, variables
            
    def graph_regulation(self, included_variables, p_vals):
        ds = self._design_space
        gma = DSDesignSpaceGMASystem(ds._swigwrapper)       
        all_variables = ds.dependent_variables + ds.independent_variables
        show_variables = [i for i in ds.dependent_variables if i not in ds.auxiliary_variables]
        show_variables += included_variables
        network_data,variable_links = self.network_data(p_vals)
        connectivity = self.connectivity()
        flux_identifiers = self.flux_identifiers()
        regulation = []
        positive_r = '[arrowhead=vee,color=gray]'
        negative_r = '[arrowhead=vee,color=gray]'
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
        
    def graph_properties(self, graph_type, size):
        properties =  'graph[ratio="fill",layout='+graph_type + ',normalize=true,size="{0},{1}"];'.format(size[0], size[1])
        properties += 'node[shape=plaintext];'
        ## properties += 'rankdir=LR;'
        return properties
        
    def subgraph_properties(self, key):
        properties = 'rank=same;'
        properties += 'color=none;'
        properties += 'edge[weight=100];'
        return properties
        
    def graph_description(self, graph_type='dot', included_variables=[],
                          p_vals=None, cmap=None, show_regulation=True,
                          size=[3.33, 2.]):
        network_data,variable_links = self.network_data(p_vals, cmap=cmap)
        regulation = self.graph_regulation(included_variables, p_vals)
        graph_string = 'digraph {'
        graph_string += self.graph_properties(graph_type, size)
        graph_string += 'subgraph {'
        for key in network_data:
            graph_string += 'subgraph cluster_'+ str(key) + ' {'
            graph_string += self.subgraph_properties(key)
            graph_string += ';'.join(network_data[key])
            graph_string += '}'
        graph_string += '}'
        if show_regulation is True:
            graph_string += 'subgraph {rankdir=LR;'#constraint=false;concentrate=true;'
            graph_string += ';'.join(regulation)
            graph_string += '}'
        graph_string += '}'
        return graph_string
    
    def graph(self, graph_type='dot', included_variables=[], 
              p_vals=None, cmap=None, show_regulation=True,
              size=[3.33,3]):
        data = {}
        graph_string = self.graph_description(graph_type=graph_type,
                                              included_variables=included_variables,
                                              p_vals=p_vals, cmap=cmap,
                                              show_regulation=show_regulation,
                                              size=size)
        data['description'] = graph_string
        if p_vals is not None:
            concentrations, [vmin, vmax] = self.flux_concentrations(p_vals)
            data['limits'] = [vmin, vmax]
        return data
             
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