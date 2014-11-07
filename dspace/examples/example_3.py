##########################
# Design Space Example 3
##########################
# This example demonstrates how the Design Space Toolbox handles
# intersecting, or overlapping, cases, which can indicate the potential
# for hysteresis.

# This example has been created from the original examples for the matlab 
# version of the design space toolbox, written by Rick Fasani.

## Import package with basic functionality
# The initial step in running the design space toolbox is importing the base
# package, and in this case the plotting functionality.
#
# Importing the plotting utilities subpackage adds plotting functions to objects
# object that are often of interest. We must also import the plotting library,
# and we make plotting interactive.

import dspace
import dspace.plotutils
from math import *
from matplotlib.pyplot import *
matplotlib.interactive(True)

## Enter Equations
# The following system has two dependent variables and two independent
# variables.  Ultimately, the design space will vary over the two
# independent variables X3 and X4.
#
#   dX1/dt = 1000*X1^2*X2^-1 + 1000*X2^-1 - X1*X3*X4
#   dX2/dt = X1^2 + 10000 - X2
#
# The equations and dependent variables are entered as strings.

f = ['X1. = alpha*X1^2*X2^-1 + alpha*X2^-1 - X1*X3*X4',
     'X2. =      X1^2 + alpha - X2']
# The next step is to construct an Equations object using this list of strings
# and specifying which variables are auxiliary and thus dependent. These are 
# specified by a list of variables names passed by the 'auxiliary_variables'
# keyword argument.

eq = dspace.Equations(f)

## Construct the System Design Space
# The previous system of equations can be analyzed using system design space.
# To construct the design space, we use an Equations object specifying our 
# system.

ds = dspace.DesignSpace(eq)

## Test for Validity
# Cases are checked to see if they are valid and the dominant terms are
# printed for each valid case.
 
valid_cases = ds.valid_cases()
cases = ds(valid_cases)
print 'Valid Cases:'
for i in cases:
    print str(i.case_number) + ':' + i.signature

## Parameter set
# We define a nominal point in parameter space.  
# The value of the individual variables at this point will not affect the plots.
# We obtain the list of parameters and independent variables from the design 
# space object.

ivar_names = ds.independent_variables

# We define a variable pool that holds symbols for all the independent variables
# and parameters of the system, by passing these names to the construction of
# the object.

pvals = dspace.VariablePool(names=ivar_names)
pvals['alpha'] = 1000
pvals['X3'] = 1  # Will be a variable on an axis, value wont be affect plot
pvals['X4'] = 1  # Will be a variable on an axis, value wont be affect plot

## Search for Intersections
# Under certain conditions, the design space can be exhaustively searched
# for all occurances of intersecting cases.  Occurances of three
# intersecting cases are especially interesting--they indicate the
# potential for bi-stability and hysteresis.  Here, the one occurance of three 
# intersecting cases (1, 2, and 4) is obtained.

intersections = ds.intersecting_cases(3, ds.valid_cases())
print intersections

##
# We make the design space object plot a 2D slice with a specified x and
# y axis and a reference parameter set. The default behavior of the plotting
# functions calculate all intersecting cases of up to 5 regions.

# Interesections of more or less cases, or specific intersections, can be draw
# by using the 'intersections' keyword argument, which receives a list of
# numbers indicating the intersecting regions that will be drawn.

fig = figure(1)
clf()
ax1 = fig.add_subplot(221)
colors = ds.draw_2D_slice(
          ax1,
          pvals,               # Pass the reference parameter set.
          'X3',                # First, indicate the x-axis variable.
          'X4',                # Second, indicate the y-axis variable.
          [1e-3, 1e3],         # Indicate the range on the x-axis.
          [1e-3, 1e3],         # Indicate the range on the y-axis.
          intersections=[1, 3] # Single cases and intersections of three.
          )
ax1.grid(True)

ax2 = fig.add_subplot(223)
colors = ds.draw_2D_slice(ax2,
                          pvals,       # Pass the reference parameter set.
                          'X3',        # First, indicate the x-axis variable.
                          'X4',        # Second, indicate the y-axis variable.
                          [1e-3, 1e3], # Indicate the range on the x-axis.
                          [1e-3, 1e3], # Indicate the range on the y-axis.
                          intersections=[3],
                          color_dict=colors
                          )
ax2.grid(True)

# The stability of the regions can be determined according to the routh
# criterion. The number of positive roots will be displayed on the z axis
# represented by a heatmap. Regions of overlap may display different values
# that correspond for different number of positive roots.  A 0 indicates
# a stable fixed point, a 1 indicates an exponentially unstable fixed point, and
# a 2 usually indicates an unstable focus for oscillations.

ax3 = fig.add_subplot(222)
cf, colors = ds.draw_2D_positive_roots(
              ax3,
              pvals,       # Pass the reference parameter set.
              'X3',        # First, indicate the x-axis variable.
              'X4',        # Second, indicate the y-axis variable.
              [1e-3, 1e3], # Indicate the range on the x-axis.
              [1e-3, 1e3], # Indicate the range on the y-axis.
              resolution=150
              )
ax3.plot([-2, 2], [0, 0], 'k.--')
ax3.grid(True)

# A 1D slice showing the number of positive roots can also be generated.  The
# number of positive roots can be displayed by different types of lines. Solid
# black lines represents zero positive roots (monostable), dotted red lines 
# represents one positive root (exponential instability) and yellow dashed 
# lines represent two positive roots (unstable focus).

ax4 = fig.add_subplot(224)
ds.draw_1D_positive_roots(ax4,
                          'log(X1)',   
                          pvals,       # Pass the reference parameter set.
                          'X3',        # First, indicate the x-axis variable.
                          [1e-2, 1e2], # Indicate the range on the x-axis.
                          )
xlabel('log(X3)')
ylabel('log(X1)')
xlim([-3, 3])
grid(True)

# We readjust the size of the axis by adding a colorbar to the plot.

ax = matplotlib.colorbar.make_axes(ax4)[0]  

# This resizes the last subplot and aligns it with the other panels in 
# the figure. We make the colorbar invisible to make the are appear empty.

ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
box(False)
show()

# The stability of individual cases can also be determined by filtering
# the cases that are drawn.  We will draw the stability of each case 
# individualy in three seperate panels. 

# We add entries to the color dictionary to make stability of the individual
# cases comparable to those of the overlapping cases. 

colors['1'] = colors['0,1']
colors['2'] = colors['0,1,2']

fig = figure(2)
clf()

ax1 = fig.add_subplot(131)
ds.draw_2D_positive_roots(ax1,
                          pvals,       # Pass the reference parameter set.
                          'X3',        # First, indicate the x-axis variable.
                          'X4',        # Second, indicate the y-axis variable.
                          [1e-3, 1e3], # Indicate the range on the x-axis.
                          [1e-3, 1e3], # Indicate the range on the y-axis.
                          resolution=150,
                          included_cases=[1],
                          color_dict=colors
                          )
ax1.set_title('Stability of Case 1:' + ds(1).signature)
ax1.grid(True)

ax2 = fig.add_subplot(132)
ds.draw_2D_positive_roots(ax2,
                          pvals,       # Pass the reference parameter set.
                          'X3',        # First, indicate the x-axis variable.
                          'X4',        # Second, indicate the y-axis variable.
                          [1e-3, 1e3], # Indicate the range on the x-axis.
                          [1e-3, 1e3], # Indicate the range on the y-axis.
                          resolution=150,
                          included_cases=[2],
                          color_dict=colors
                          )
ax2.set_title('Stability of Case 2:' + ds(2).signature)
ax2.grid(True)

ax3 = fig.add_subplot(133)
ds.draw_2D_positive_roots(ax3,
                          pvals,       # Pass the reference parameter set.
                          'X3',        # First, indicate the x-axis variable.
                          'X4',        # Second, indicate the y-axis variable.
                          [1e-3, 1e3], # Indicate the range on the x-axis.
                          [1e-3, 1e3], # Indicate the range on the y-axis.
                          resolution=150,
                          included_cases=[3],
                          color_dict=colors
                          )
ax3.set_title('Stability of Case 3:' + ds(3).signature)
ax3.grid(True)

