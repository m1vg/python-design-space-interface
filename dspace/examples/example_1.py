##########################
# Design Space Example 1
##########################
#
# This simple example demonstrates the basic steps and functions necessary
# to construct and analyze a design space using the Design Space Toolbox V2
# python package.
#
# This example has been created using the original examples for the matlab 
# version of the design space toolbox, written by Rick Fasani.

## Import package with basic functionality
# The initial step in running the design space toolbox is importing the base
# package. Importing this package adds all the objects and functions to analyze
# a system WITHOUT importing plotting utilities. 

import dspace

## Enter Equations
# The following system has two dependent variables and two independent
# variables and one parameter, alpha.  Ultimately, the design space will vary 
# over the two independent variables X3 and X4 with a specified value of alpha.
#
# dX1/dt =     alpha*X1*X2*X3*X4 + X1*X2 - X1
# dX2/dt = alpha^-1*X1*X2*X3*X4 + X1*X2 - X2
#
# Each equations in string format is entered forming a list of equations, where
# the left hand side of each equation is indicated by either the 
# time-differentiated variable using dot notation for differentiation, or zero
# for constraint equations.  All kinetic orders must not be variables (although
# symbols may be replaced by values when constructing a system design space
# object).

f = ['X1. =     alpha*X1*X2*X3*X4 + X1*X2 - X1',
     'X2. = alpha^-1*X1*X2*X3*X4 + X1*X2 - X2']
 
# The next step is to construct an Equations object using this list of strings.

eq = dspace.Equations(f)

## Construct the System Design Space
# The previous system of equations can be analyzed using system design space.
# To construct the design space, we use an Equations object specifying our 
# system.

ds = dspace.DesignSpace(eq)

##
# The term signature, or number of positive and negative terms in each
# equation, is important throughout the analysis.

print('Terms: ' + ds.signature)

## Retrieve Cases
# A design space is completely characterized by a set of dominant S-systems
# and the conditions for their dominance.  
#
# A case can be retrieved by the system design space object by calling it with
# a case number of interest.

case_1 = ds(1)

# Alternatively a list of case can be retrieved, and the list can be unpacked.

cases = ds([2, 3, 4])
print(cases[0].equations)

case_2, case_3, case_4 = cases

# Cases can be retrieved using their case signature.

case_1121 = ds('1121', by_signature=True)

##
# The dominant S-systems and conditions for each case can be printed to the
# command window.  Here, only the first case is shown.

print('Equations: '+str(case_1.equations))
print('Conditions: '+str(case_1.conditions_log))
print('Boundaries: '+str(case_1.boundaries_log))

## Test for Validity
# Some sets of boundaries are mutually exclusive, or invalid.  Under certain
# conditions, the validity of boundaries can be checked automatically.
#
# We obtain a list of case numbers that are valid.
 
valid_cases = ds.valid_cases()
print('All valid cases: ' + str(valid_cases))

# In some cases we want to restrict parameter values to a range or a fixed value
# and test for validity within these constraints.

valid_cases = ds.valid_cases(p_bounds={'alpha':10,
                                       'X3':[1e-3, 1e3],
                                       'X4':[1e-3, 1e3]
                                       })
print('All valid cases within bounds: ' + str(valid_cases))

##
# The dominant terms in each case can be printed to the command window.
# Here, only the valid cases are shown.

# We loop over all the valid cases by their case number
for case_number in valid_cases:
    # Get the case object
    case = ds(case_number)
    # Print the case signature
    print(str(case.case_number) + ':' + case.signature)

## Plot

## Import package with plotting functionality
# The basic functionality does not include plotting utilities. However, included
# in the dspace package are plotting utilities that rely on NumPy, SciPy and
# Matplotlib.
# 
# Importing the plotting utilities subpackage adds plotting functions to objects
# object that are often of interest. We must also import the plotting library,
# and we make plotting interactive.

import dspace.plotutils
from matplotlib.pyplot import *
matplotlib.interactive(True)

# Regions in design space are made up of one or more overlapping cases.  To
# visualize the regions, individual cases are checked to see if they lie
# inside the boundaries of the plotting and if so, their vertices are enumerated
# and polytopes are drawn.
#
# First, the a nominal point in parameter space must be designated.  
# The value of the individual variables at this point will not affect the plots.
# We obtain the list of parameters and independent variables from the design 
# space object.

ivar_names = ds.independent_variables

# We define a variable pool that holds symbols for all the independent variables
# and parameters of the system, by passing these names to the construction of
# the object.

pvals = dspace.VariablePool(names=ivar_names)
pvals['alpha'] = 10
pvals['X3'] = 1  # Will be a variable on an axis, value wont be affect plot
pvals['X4'] = 1  # Will be a variable on an axis, value wont be affect plot

##
# We make the design space object plot a 2D slice with a specified x and
# y axis and a reference parameter set.


clf()
fig = gcf()
ax = gca()
colors = ds.draw_2D_slice(ax,
                          pvals,       # Pass the reference parameter set.
                          'X3',        # The x-axis variable.
                          'X4',        # The y-axis variable.
                          [1e-3, 1e3], # Range on the x-axis.
                          [1e-3, 1e3], # Indicate the range on the y-axis.
                          )
grid(True)
show()

##
# We make the design space object plot an interactive 2D slice with a specified 
# x and y axis, a reference parameter set and a dictionary with parameter, range
# pairs indicating sliders to exlore parameter space.

colors['3'] = (1, 1, 0)

ds.draw_2D_slice_interactive(
          pvals,                          # Pass the reference parameter set.
          'X3',                           # The x-axis variable.
          'X4',                           # The y-axis variable.
          [1e-3, 1e3],                    # The range on the x-axis.
          [1e-3, 1e3],                    # The range on the y-axis.
          {'alpha':[1e-5, 1e5]},          # Specify slider parameters and ranges
          color_dict = colors
          )
