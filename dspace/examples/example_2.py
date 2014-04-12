##########################
# Design Space Example 2
##########################
# This example is similar to the previous example, but demonstrates more
# advanced functions of the Design Space Toolbox.

# This example has been created from the original examples for the matlab 
# version of the design space toolbox, written by Rick Fasani.

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
# The following system is similar to the previous example, but with
# multiple negative terms.  It still has two dependent variables and two
# independent variables, and the design space will vary over the two
# independent variables X3 and X4.
#
# dX1/dt =     10*X1*X2*X3*X4 + X1*X2 -     10*X1*X3*X4^-1 - X1
# dX1/dt = (1/10)*X1*X2*X3*X4 + X1*X2 - (1/10)*X2*X3*X4^-1 - X2
#
# Each equations in string format is entered forming a list of equations, where
# the left hand side of each equation is indicated by either the 
# time-differentiated variable using dot notation for differentiation, or zero
# for constraint equations.  An example of a algebraic constraint is shown for
# the system above.
#
# dX1/dt =     10*X5 + X1*X2 -     10*X1*X3*X4^-1 - X1
# dX1/dt = (1/10)*X5 + X1*X2 - (1/10)*X2*X3*X4^-1 - X2
#      X5 = X1*X2*X3*X4
#

f = ['X1. =     alpha*X5 + X1*X2 -     alpha*X1*X3*X4^-1 - X1',
     'X2. = (1/alpha)*X5 + X1*X2 - (1/alpha)*X2*X3*X4^-1 - X2',
     '  X5 = X1*X2*X3*X4']
 
# The next step is to construct an Equations object using this list of strings
# and specifying which variables are auxiliary and thus dependent. These are 
# specified by a list of variables names passed by the 'auxiliary_variables'
# keyword argument.

eq = dspace.Equations(f, auxiliary_variables=['X5'])

## Construct the System Design Space
# The previous system of equations can be analyzed using system design space.
# To construct the design space, we use an Equations object specifying our 
# system.

ds = dspace.DesignSpace(eq)

## Test for Validity
# We choose X3 and X4 as are axes for a 2D plot, and must specify a value for
# alpha. We obtain the list of valid cases with a slice through design space
# with alpha = 10, 1e-3 < X3 < 1e3 and 1e-3 < X4 < 1e3. 
 
valid_cases = ds.valid_cases(p_bounds={'alpha':10, 
                                       'X3':[1e-3, 1e3],
                                       'X4':[1e-3, 1e3]})

print 'All valid cases in slice: ' + str(valid_cases)

## Plot
# We define a nominal point in parameter space.  
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

## Global color scheme for use across multiple slices.
# In this example we will plot many different slices of the same dsign space.
# To get consistent colors for the different regions, we will use a global color
# scheme with all the valid cases included.

all_valid_cases = ds.valid_cases()
colors = dict()
for i in range(0, len(all_valid_cases)):
    colors[str(all_valid_cases[i])] = matplotlib.cm.hsv(float(i)/len(all_valid_cases))

temp = colors['6']
colors['6'] = colors['11']
colors['11'] = temp
##
# We make the design space object plot a 2D slice with a specified x and
# y axis and a reference parameter set. To make the colors consistent across 
# the figures we have included the colors from the previous plot by passing a 
# color dictionary to the 'color_dict' keyword argument in the drawing method.

fig = figure(1)
clf()
ax = gca()
ds.draw_2D_slice(ax,
                 pvals,             # Pass the reference parameter set.
                 'X3',              # First, indicate the x-axis variable.
                 'X4',              # Second, indicate the y-axis variable.
                 [1e-3, 1e3],       # The range on the x-axis.
                 [1e-3, 1e3],       # The range on the y-axis.
                 color_dict=colors, # Specify a color dictionary for the cases.
                 )
sca(ax)

## Obtaining the vertices in 2D
# The verticesfor the 2D slice are obtained for case 4, the diamond in the
# center.

case4 = ds(4)
vertices = case4.vertices_2D_slice(pvals, 'X3', 'X4')
print vertices

## Create a Bounding Box in higher domensions
# A bounding box is calculated for case 4. The result is given as lower and 
# upper bounds, respectively, on the parameters, X3 and X4 (in Cartesian 
# coordinates).

box = case4.bounding_box(log_out=True)
print box

# A bounding box can be obtained with specific constraints on the parameters,
# by passing parameter values or parameter ranges.


box = case4.bounding_box(p_bounds={'alpha':10, 
                                   'X3':[1e-3, 1e3],
                                   'X4':[1e-3, 1e3]}, 
                         log_out=True);
print box

##
# The box is plotted in the previous figure, where the bounds (in
# logarithmic coordinates) can be visually confirmed.


x_vertices = [min(box['X3']), min(box['X3']),max(box['X3']), max(box['X3'])]
y_vertices = [min(box['X4']), max(box['X4']), max(box['X4']), min(box['X4'])]
fill(x_vertices, y_vertices,facecolor='none', edgecolor='k', linewidth=2.)

## Measure Tolerance
# An operating point is set such that it is in case 4.

pvals['X3'] = 2
pvals['X4'] = 1

# The tolerance from this
# point to the boundaries of the region is measured for each orthogonal 
# direction.  The result is given as fold change down and fold change up,
# respecively (in Cartesian coordinates).

tol = case4.measure_tolerance(pvals, log_out=True);
print tol

# The operating point is plotted and the tolerance (in logarithmic coordinates) 
# can be visually confirmed.  We make a new plot with only the case of interest,
# and overlay the operating point and measured tolerances.

figure(2)
clf()
ax = gca()
case4.draw_2D_slice(ax, pvals, 'X3', 'X4', [1e-3, 1e3], [1e-3, 1e3], 
                    facecolor=colors['4'], edgecolor='none')
xlim(-3, 3)
ylim(-3, 3)


# Plot the operating point as a black point.

plot(log10(pvals['X3']), log10(pvals['X4']), 'k.')

# Plot the tolerances as two lines using list enumeration.
plot([log10(pvals['X3'])]*2,
     [i+log10(pvals['X4']) for i in tol['X4']], 'k')
plot([i+log10(pvals['X3']) for i in tol['X3']],
     [log10(pvals['X4'])]*2, 'k')

## Plot Multiple Cases
# The design space is plotted for two cases: alpha = 10 and alpha = 1/10.

fig = figure(3)
clf()

ax1 = fig.add_subplot(121)
pvals['alpha'] = 10
ds.draw_2D_slice(ax1, pvals, 'X3', 'X4', [1e-3, 1e3], [1e-3, 1e3],
                 color_dict=colors)
sca(ax1)
xlabel(r'$\log_{10}(X_3)$');
ylabel(r'$\log_{10}(X_4)$');
title(r'$\alpha = 10$');
grid(True);

ax2 = fig.add_subplot(122)
pvals['alpha'] = 1./10
ds.draw_2D_slice(ax2, pvals, 'X3', 'X4', [1e-3, 1e3], [1e-3, 1e3], 
                 color_dict=colors)
sca(ax2)
xlabel(r'$\log_{10}(X_3)$');
ylabel(r'$\log_{10}(X_4)$');
title(r'$\alpha = 1/10$');
grid(True);

## Plot 3D Slice 
# A plot with alpha changing from 10 to 1/10, X3, and X4 produce a 3-D image 
# that can be constructed and manipulated using th plot utilities. First we must
# import the 3D plotting component of matplotlib.


from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as a3

##
# We make the design space object plot the 3D slice a specified x, y and z axes
# with reference parameter set.

fig = figure(4)
clf()

ax = fig.add_subplot('111', projection='3d')       # We define the plot axes where it will be drawn.
ds.draw_3D_slice(ax,
                 pvals,       # Pass the reference parameter set.
                 'X3',        # First, indicate the x-axis variable.
                 'X4',        # Second, indicate the y-axis variable.
                 'alpha',     # Thir, indicate the z-axis variable.
                 [1e-3, 1e3], #Indicate the range on the x-axis.
                 [1e-3, 1e3], #Indicate the range on the y-axis.
                 [1e-1, 1e1], #Indicate the range on the z-axis.
                 color_dict=colors
                )
