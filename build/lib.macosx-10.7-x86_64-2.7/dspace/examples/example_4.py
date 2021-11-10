##########################
# Design Space Example 4
##########################
# This simple example demonstrates the analysis of the dominant S-system
# in each region using the Design Space Toolbox.

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
# dX1/dt =     10*X1*X2*X3*X4 + X1*X2 - X1
# dX1/dt = (1/10)*X1*X2*X3*X4 + X1*X2 - X2
#
# Here, the system is input, parsed, and analyzed.
# The equations and dependent variables are entered as strings.

f = ['X1. =     10*X1*X2*X3*X4 + X1*X2 - X1',
     'X2. = (1/10)*X1*X2*X3*X4 + X1*X2 - X2']

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
pvals['X3'] = 1  # Will be a variable on an axis, value wont be affect plot
pvals['X4'] = 1  # Will be a variable on an axis, value wont be affect plot

## Plot Region
# Here, the regions are plotted as in the other examples. 
# We make the design space object plot a 2D slice with a specified x and
# y axis and a reference parameter set.

fig = figure(1)
clf()
ax1 = gca()
colors = ds.draw_2D_slice(
          ax1,
          pvals,               # Pass the reference parameter set.
          'X3',                # First, indicate the x-axis variable.
          'X4',                # Second, indicate the y-axis variable.
          [1e-3, 1e3],         # Indicate the range on the x-axis.
          [1e-3, 1e3],         # Indicate the range on the y-axis.
          )
ax1.grid(True)
show()

## Sample the Steady States
# In every region, the dominant S-system can be analyzed using well-known
# methods from Biochemical Systems Theory and the Power Law Formalism.
# Here, the steady state for X1 and X2 in log space is calculated at each point
# as defined by the resolution of the plot.

fig = figure(2)
clf()

ax1 = fig.add_subplot(221)
colors = ds.draw_2D_ss_function(
          ax1,
          'log(X1)',           # Single-value function to plot as a string.
          pvals,               # Pass the reference parameter set.
          'X3',                # First, indicate the x-axis variable.
          'X4',                # Second, indicate the y-axis variable.
          [1e-3, 1e3],         # Indicate the range on the x-axis.
          [1e-3, 1e3],         # Indicate the range on the y-axis.
          resolution=100,      # Approximately a 100x100 plot.
          )
ax1.plot([-2, 2], [0, 0], 'k.--')
ax1.grid(True)

ax2 = fig.add_subplot(222)
colors = ds.draw_2D_ss_function(
          ax2,
          'log(X2)',           # Single-value function to plot as a string.
          pvals,               # Pass the reference parameter set.
          'X3',                # First, indicate the x-axis variable.
          'X4',                # Second, indicate the y-axis variable.
          [1e-3, 1e3],         # Indicate the range on the x-axis.
          [1e-3, 1e3],         # Indicate the range on the y-axis.
          resolution=100,      # Approximately a 100x100 plot.
          )
ax2.plot([-2, 2], [0, 0], 'k.--')
ax2.grid(True)
show()

# A 1D slice showing the concentration on the y-axis can also be drawn.

ax3 = fig.add_subplot(223)
ds.draw_1D_ss_function(ax3,
                       'log(X1)',   
                       pvals,       # Pass the reference parameter set.
                       'X3',        # First, indicate the x-axis variable.
                       [1e-2, 1e2], # Indicate the range on the x-axis.
                       )
xlabel('log(X3)')
ylabel('log(X1)')
xlim([-2, 2])
grid(True)

ax4 = fig.add_subplot(224)
ds.draw_1D_ss_function(ax4,
                      'log(X2)',   
                      pvals,       # Pass the reference parameter set.
                      'X3',        # First, indicate the x-axis variable.
                      [1e-2, 1e2], # Indicate the range on the x-axis.
                      )
xlabel('log(X3)')
ylabel('log(X2)')
xlim([-2, 2])
grid(True)

# We add a new axis to subplot 2 and 3 to resize them and align them with
# the other figures.

ax = matplotlib.colorbar.make_axes(ax3)[0]  
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
box(False)
show()

ax = matplotlib.colorbar.make_axes(ax4)[0]  
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
box(False)
show()

# Individual cases can be drawn individually by performing the same operation
# drawing method of the case object.  The z-limits of the case will not span
# the complete range of the original plot, thus we apply the use the same 
# limits for comparison.

case2 = ds(2)

fig = figure(3)
clf()

ax1 = gca()
colors = case2.draw_2D_ss_function(
          ax1,
          'log(X2)',           # Single-value function to plot as a string.
          pvals,               # Pass the reference parameter set.
          'X3',                # First, indicate the x-axis variable.
          'X4',                # Second, indicate the y-axis variable.
          [1e-3, 1e3],         # Indicate the range on the x-axis.
          [1e-3, 1e3],         # Indicate the range on the y-axis.
          resolution=100,      # Approximately a 100x100 plot.
          zlim=[-5, 0],
          )
ax1.grid(True)

