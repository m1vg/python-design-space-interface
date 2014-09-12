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
matplotlib.interactive(True)

## Enter Equations
# The following system has three dependent variables and two independent
# variables.  Ultimately, the design space will vary over the two
# independent variables X3 and X4 and two parameters, alpha and beta.
#
# dX1/dt =          alpha*X5 + X4 - X1
# dX2/dt = X1 - X2 - X2*X6 - beta*X2^2
# dX3/dt =                     X2 - X3
# dX4/dt =             X3 + X2*X6 - X4
#
# Here, the system is input, parsed, and analyzed.
# The equations and dependent variables are entered as strings.

f = ['X1. = alpha*X5 + 2*k41*X4 - 2*X1^2', 
 'X2. = X1^2 + X3*X7 - k23*X2 - X2*X6',# - beta*X2^2
 'X3. = k23*X2 - k34*X3 - X3*X7',
 'X4. = k34*X3 + X2*X6 - k41*X4 - beta*X4^2']



## f = ['X1. = alpha*X5 + 2*X4 - 2*X1^2', 
##  'X2. = X1^2 - X2 - X2*X6 - beta*X2^2',
##  'X3. = X2 - X3',
##  'X4. = X3 + X2*X6 - X4']

## f = ['X1. = alpha*X5 + 2*k41*X4 - 2*X1^2', 
##  'X2. = X1^2 - k23*X2 - X2*X6 - beta*X2^2',
##  'X3. = k23*X2 - k34*X3',
##  'X4. = k34*X3 + X2*X6 - k41*X4']

# The next step is to construct an Equations object using this list of strings.

eq = dspace.Equations(f)

## Construct the System Design Space - with and without resolving cycles.
# The previous system of equations can be analyzed using system design space.
# First, we will define a system design space in the way we have defined it
# in the previous examples, without resolving any cycles.

ds_nocycles = dspace.DesignSpace(eq, name='Without Cycles')

# Then, we will resolve the cycles in a second design space object and we will 
# compare it with the design space where the cycles are not resolved. To resolve
# cycles, the design space is constructed with the 'resolve_cycles=true' 
# keyword argument.

ds_cycles = dspace.DesignSpace(eq, resolve_cycles=True, name='With Cycles')#,constraints=constraints)

## Test for Validity
# We calculate the valid cases without and with resolving cycles. As we can see
# the valid cases are not the same for the two methods.  The set of valid cases
# for the design space without resolving cycles is a subset of those found when
# resolving cycles. 

valid_no_cycles = ds_nocycles.valid_cases()
valid_cycles = ds_cycles.valid_cases()
## 
## print 'All valid cases without resolving cycles: ' + str(valid_no_cycles)
## print 'All valid cases with resolving cycles: ' + str(valid_cycles)
## 
## print 'Is subset: ' + str(set(valid_no_cycles).issubset(set(valid_cycles)))
## 
## ## Cyclical cases
## # The additional cases in the design space with resolving cycles are cases
## # where the dominant terms create an underdetermined system that does not have
## # a solution. These singular cases are resolved by including subdominant 
## # processes. For each cyclical case there may be a number of dominant subcases.
## 
## # We first determine the cases that are cyclical by iterating over all the valid
## # cases.
## 
## for case in ds_cycles(valid_cycles):
##     print str(case.case_number)+':'+case.signature+ ' cyclical? ' + str(case.is_cyclical)
## 
## # A cyclical case is similar to a design space in that it has subcases. Of these
## # subcases not all of them are valid and we can determine the subcases that are
## # valid.
## 
## valid_subcases = ['7_'+ str(i) for i in ds_cycles(7).valid_subcases()]
## valid_subcases += ['10_'+ str(i) for i in ds_cycles(10).valid_subcases()]
## 
## print valid_subcases
## 
## # A subcase can be obtained from the case object similar to how a case is 
## # obtained from a design space.
## 
## case10 = ds_cycles(10)
## case10_6 =  case10(6)
## 
## # Subcases can be obtained directly from the design space object, by passing 
## # a string with the following format: 
## #    <case number>_<subcase number>
## 
## case10_6 = ds_cycles('10_6')
## 
## # A subcase is an instance of the Case class and thus all methods that are
## # available for a case is available to a subcase. As an example, we will obtain
## # a parameter set that makes Case 7 Subcase 6 the dominant case.
## 
## pvals = case10_6.valid_parameter_set()
## print pvals
## 
## ## Cyclical subcases
## # When resolving a cyclical case with subdominant fluxes, it is possible to 
## # create a new cycle that again must be resolved. The current implementation 
## # does not solve all possible subcycles [under development], but can solve some.
## # This is the case for Case 7 subcase 3. We inspect the equations of the subcase
## # to identify the cycle directly.
## 
## case10_3 = ds_cycles('10_3')
## 
## print case10_3.name + ' is cyclical? ' + str(case10_3.is_cyclical)
## for equation in case10_3.equations:
##     print equation
## 
## # There is a new cycle between the equation for X2 and the equation for X3.
## # we can print the valid subcases for case 10 subcase 3 much in the same way
## # we found and printed the valid subcases for case 10.
## 
## print  ['10_3_'+i for i in case10_3.valid_subcases()]
## 
## # A more convenient way to obtain all the subcases, and their subcases, is to
## # pass a list of cases and remove all the cyclical cases and replacing them 
## # with their respective subcases.
## 
## print ds_cycles.cycles_to_subcases(valid_cycles)
## 
## ## All valid cases by resolving cycles.
## # The list of cases obtained by expanding the cycles includes both valid
## # and invalid subcases. It is more often useful to restrict the list of 
## # cases and subcases to those that are valid. A list of all valid cases and
## # expanding this list to all valid subcases is obtained using the 'valid_cases'
## # method for the design space object.
## 
## print ds_cycles.valid_cases(expand_cycles=True)
## 
## # The list of cases and subcases bounded in parameter space may also be
## # obtained.
## 
## print ds_cycles.valid_cases(expand_cycles=True, p_bounds={'alpha':100,
##                                                            'beta':1e-1,
##                                                            'X5':1e-1,
##                                                            'X6':[1e-1, 1e0]})

## Plot
# We define a nominal point in parameter space and plot the design spaces 
# without and with resolving the cycles.

pvals=dspace.VariablePool(names=ds_cycles.independent_variables)
pvals['alpha'] = 10
pvals['beta'] = 0.1
pvals['X7'] = 1.5


ds_cycles.draw_2D_slice_interactive(pvals,             # Pass the reference parameter set.
                                    'X5',              # First, indicate the x-axis variable.
                                    'X6',              # Second, indicate the y-axis variable.
                                    [1e-3, 1e3],       # The range on the x-axis.
                                    [1e-3, 1e3],       # The range on the y-axis.
                                    {'X7':[1e-5, 1e5]}
                                    )

fig = figure(1)
clf()

# We plot the design space with cycles first to get the colors of all the
# cases and subcases. These colors are used to make the plot without resolving
# the cycles to match the colorbars.

ax = fig.add_subplot('122')       # We define the plot axes.
colors = ds_cycles.draw_2D_slice(ax,
                             pvals,             # Pass the reference parameter set.
                             'X5',              # First, indicate the x-axis variable.
                             'X6',              # Second, indicate the y-axis variable.
                             [1e-3, 1e3],       # The range on the x-axis.
                             [1e-3, 1e3],       # The range on the y-axis.
                             )
ax.set_title(ds_cycles.name)

ax = fig.add_subplot('121')       # We define the plot axes.
ds_nocycles.draw_2D_slice(ax,
                      pvals,             # Pass the reference parameter set.
                      'X5',              # First, indicate the x-axis variable.
                      'X6',              # Second, indicate the y-axis variable.
                      [1e-3, 1e3],       # The range on the x-axis.
                      [1e-3, 1e3],       # The range on the y-axis.
                      color_dict=colors
                      )
ax.set_title(ds_nocycles.name)


## 
## for i in ['X1', 'X2', 'X3', 'X4']:
##     fig=figure()
##     ds_cycles.draw_2D_ss_function(gca(),
##                                   'log(V_'+i+')',
##                                   pvals, 
##                                   'X5',              # First, indicate the x-axis variable.
##                                   'X6',              # Second, indicate the y-axis variable.
##                                   [1e-3, 1e3],       # The range on the x-axis.
##                                   [1e-3, 1e3],       # The range on the y-axis.
##                                   )
## 


