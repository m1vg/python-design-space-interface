##########################
# Design Space Example 6
##########################

# WARNING!! This example does not currently work with the develop version of
# the C toolbox. Until further notice, running this file will cause errors to
# appear.

# This example is similar to previous examples, but demonstrates the analysis
# of systems with cycles. A cycle that is dominant creates singularity in 
# conventional design space analysis.  These singularities are resolved 
# automatically by assuming subdominant fluxes.

## Import package with basic functionality
# The initial step in running the design space toolbox is importing the base
# package, and in this case the plotting functionality.
#
# Importing the plotting utilities subpackage adds plotting functions to objects
# object that are often of interest. We must also import the plotting library,
# and w make plotting interactive.
import dspace
import dspace.plotutils
from math import *
from matplotlib.pyplot import *


f = ['X1. = a1 - X1',
     'X2. = X1 - X2',
     'X3. = X2 - X3*X4^-1',
     'X4 = X3 + b1',
     ]

f = ['X1. = X1^3*X2^2 - 6*X1^2*X2 - b1*X1',
     'X2. = 11*X2 - 6*X1^-1']

fe = ['X1. = a1*X1^3*X2^2 - 6*X1^2*X2 + epsilon^5*X1^-5  - epsilon^5*X1^5',
     'X2. = 11*X2 - b1*6*X1^-1 + epsilon^5*X2^-5 - epsilon^5*X2^5']

f = ['X1. = a1*X1^3*X2^2 - 6*X1^2*X2 - epsilon*X1',
     'X2. = 11*X2 - b1*6*X1^-1']


eq = dspace.Equations(f)

ds = dspace.DesignSpace(eq)

eq = dspace.Equations(fe)

dse = dspace.DesignSpace(eq, constraints=['epsilon < 0.000000001'])

pvals=dspace.VariablePool(names=ds.independent_variables)

pvals_e=dspace.VariablePool(names=dse.independent_variables)
pvals['epsilon']=1e-20
pvals_e['epsilon']=1e-20
## pvals_e['epsilon2']=1e-10

ics = {i:1e-10 for i in ds.dependent_variables}

cdict={'1': (1.0, 0.0, 0.16, 1.0),
       '3': (0.36036036036036034, 1.0, 0.0, 1.0),
       '4': (0.0, 0.56159420289855033, 1.0, 1.0)
       }

ax=subplot(1,2,1)
for i in xrange(1, ds.number_of_cases+1):
    c=ds(i)    
    c.draw_2D_phase_portrait(gca(), 
                             ics, 
                             pvals, 'X1', 'X2', [1e-20, 1e20], [1e-20, 1e20], 
                             resolution=20, show_designspaces=True, 
                             color_dict=cdict)

ax=subplot(1,2,2)
for i in xrange(1, dse.number_of_cases+1):
    c=dse(i)   
    c.draw_2D_phase_portrait(gca(), 
                             ics, 
                             pvals_e, 'X1', 'X2', [1e-10, 1e15], [1e-10, 1e15], 
                             resolution=20, show_designspaces=True, 
                             color_dict=cdict)

figure()
for i in [1, 2]:
    c=ds(i)
    c.draw_2D_phase_portrait(gca(), 
                             ics, 
                             pvals, 'X1', 'X2', [1e0, 10**0.5], [1e0, 10**0.5], 
                             resolution=20, show_designspaces=True, 
                             color_dict=cdict)