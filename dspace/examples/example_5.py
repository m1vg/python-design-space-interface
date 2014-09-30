##########################
# Design Space Example 5
##########################
# This simple example demonstrates the use of the runutils module for routine
# analysis of design space objects.

## Import the runutils module of the dspace package.
# The initial step in using the runutils module is importing its functionality
# into the namespace.  From this module we only use to Input class that 
# holds the input criteria for a particular design space.

from dspace.runutils import Input

## Specifying the input criteria.
# The Input object must accept a list of equations that specifies the
# model.  Several other keyword arguments are allowed, as described by the 
# by the Input class help documentation.

help(Input)

## Define the system and analysis criteria.
# The following system has two dependent variables and one independent
# variable.  Ultimately, the design space will vary over the
# independent variable X3 and the parameter a1.
#
# dX1/dt =      a1 + a2*X3 - b1*X1
# dX2/dt = b1*X1 + a2*X3^2 - b2*X2
#
# Next, we show an example for a 2 variable system.  This example shows how 
# most of the plotting and printing utilities are incorporated as
# arguments to the Input class.

Input(['X1. = a1 + a2*X3 - b1*X1',
       'X2. = a3*X1 + a2*X3^2 - b2*X2'], # Specify the system of equations.
      auxiliary_variables=[],            # Identify auxiliary variables.
      xaxis='X3',                        # Specify the x-axis.
      yaxis='a1',                        # Specify the y-axis.
      parameters={'a1':1, 
                  'a2':1, 
                  'b1':1e-3, 
                  'b2':1e-3,
                  'a3':1e-3,
                  'X3':100
                   },                    # Input the nominal parameter set.
      print_valid_cases=True,            # Specify if will print valid cases.
      x_range=[1e-3, 1e3],               # Range of the x-axis.
      y_range=[1e-3, 1e3],               # Range of the y-axis.
      centered_axes=True,                # Axis are centered about parameters.
      # Include interactive plot with 4 sliders: one for each a2, b2, b1 and a3.
      plot_interactive={'a2':[1e-5, 1e5],   
                        'b2':[1e-10, 1e10],
                        'b1':[1e-10, 1e10], 
                        'a3':[1e-5, 1e5]},
      
      plot_designspace=True,             # Should draw design space plot.
      plot_steady_states=['X1', 'X2'],   # Concentrations to plot.
      plot_fluxes=['X1'],                # Fluxes to plot.
      plot_stability=True,               # Should draw stability plot.
      plot_functions=['log(V_X1/X2)'],   # Arbitrary functions to plot.
      intersections=[1, 3]
      )  

## Specifying a subset of cases to draw.
# The same system can be analyzed, except a subset of cases can be selected for 
# drawing. We can specify the cases that will be drawn by using the `draw_cases`
# keyword argument.
Input(['X1. = a1 + a2*X3 - b1*X1',
       'X2. = b1*X1 + a2*X3^2 - b2*X2'],
      auxiliary_variables=[],
      xaxis='X3',
      yaxis='a2',
      #get_parameters=4,           # Automatically gets a parameter set for case 4.
      get_parameters=':2121',      # Automatically gets a parameter set for case with signature :2121. 
      print_valid_cases=True,
      x_range=[1e-3, 1e3],
      y_range=[1e-3, 1e3],
      centered_axes=True,
      plot_designspace=True,
      plot_steady_states=['X1'],
      plot_fluxes=['X1'],
      plot_stability=True,
      plot_functions=['X1'],
      draw_cases=[':11*1', 3],    # Can use signature (w/ wildcards) to get specific sets of cases
      zlim=[-10, 10] # Explicitly sets z-lim of steady state, flux and function plots
      )

