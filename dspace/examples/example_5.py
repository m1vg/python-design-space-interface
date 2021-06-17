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

Input(['X1. = a1^n + a2*X3 - b1*X1',
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
      parameter_dict={'n':1},
      plot_designspace=True,             # Should draw design space plot.
      plot_steady_states=['X1'],         # Concentrations to plot.
      plot_log_gains=[('X1','X3'),
                      ('X2','X3')],      # Log gain/sensitivity plot.
      plot_fluxes=['X1'],                # Fluxes to plot.
      plot_stability=True,               # Should draw stability plot.
      plot_functions=['log(V_X1/X2)',    # Arbitrary functions to plot.
                      '$L_X1_X3'],       # Log gains by using this notation: $L_<dependent>_<independent>.
      intersections=[1, 3],
      show_vertices=[2],                 # Show vertices for case 2
      vertex_font_size=12,               # Set the font size for vertices
      latex_symbols={'a1':r'\alpha_1',   
                     'a2':r'\alpha_2',
                     'a3':r'\alpha_3',
                     'b1':r'\beta_1',
                     'b2':r'\beta_2', 
                     'X1':r'X_1',
                     'X2':r'X_2',
                     'X3':r'X_3'},       # latex representation for parameters and variables
      resolution=100,
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
      #get_parameters=':1121',     # Automatically gets a parameter set for case with signature :2121. 
      #get_parameters=2,           # Automatically gets a parameter set for case number 2. 
      get_parameters=[2, 3],       # Automatically gets a parameter set for cases 2 and 3 [intersection or co-localization].
      colocalize_cases=True,       # Indicates get_parameters for multiple cases is case co-localization or intersection.
      #minimize_function='X1',     # Determined parameter set minimizes this objective function.
      maximize_function='X4',      # Determined parameter set minimizes this objective function.
      objective_bounds={'a1':0.1, 
                        'a2':[1e-3, 1e0],
                        'a3':0.001,
                        'X3':[1e-4, 1e4], 
                        'b1':0.45,
                        'b2':0.45
                        },
      print_valid_cases=True,
      x_range=[1e-4, 1e4],
      y_range=[1e-3, 1e0],
      centered_axes=False,
      plot_designspace=True,
      plot_steady_states=['X2'],
 
      ## plot_fluxes=['X2'],
      ## plot_stability=True,
      ## plot_functions=['X1'],
      draw_cases=[':11*1', 3],    # Can use signature (w/ wildcards) to get specific sets of cases
      zlim=[-1.5, 1.5] # Explicitly sets z-lim of steady state, flux and function plots
      )
