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
matplotlib.interactive(False)

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

system_sets = [['X1. = alpha*X5 + 2*k41*X4 + 4*gamma*X2^2 - 2*X1^2', 
                'X2. = X1^2 + k32*X3 - X2*X8 - X2*X7 - 2*gamma*X2^2',
                'X3. = X2*X8 - k34*X3*X6 - k32*X3',
                'X4. = k34*X3*X6 + X2*X7 - k41*X4 - 2*beta*X4^2',
                'X6 = 1 + X3',
                ],
               ## ['X1. = X5 + 2*k21*X2 - 2*k12*X1^2 - b1*X1',
               ##  'X2. = X6 + k12*X1^2 + k32*X3 - k21*X2 - k23*X2 - b2*X2',
               ##  'X3. = X7 + k23*X2 + k43*X4^2 - k32*X3 - k34*X3 - b3*X3',
               ##  'X4. = X8 + 2*k34*X3 - 2*k43*X4^2 - b4*X4'
               ##  ],
               ## ['X1. = 2*k21*X2 - 2*k12*X1^2*X3^-1 - b1*X1',
               ##  'X2. = X6 + k12*X1^2*X3^-1 + k32*X3 + k42*X4 - k21*X2 - k23*X2 - k24*X2',
               ##  'X3. = k23*X2 + k43*X4 - k32*X3 - k34*X3',
               ##  'X4. = X7 + k24*X2 + k34*X3 + k54*X5*X3^-1 - k42*X4 - k43*X4 - k45*X4',
               ##  'X5. = k45*X4 - k54*X5*X3^-1 - b5*X5'
               ##  ],
               ## ['X1. = a1 + 2*k21*X2 - 2*k12*X1^2',
               ##  'X2. = k12*X1^2 + X3*X7 + X4*X6 - k21*X2 - X2*X6 - X2*X7',
               ##  'X3. = X2*X6 + X4*X7 - X3*X6 - X3*X7',
               ##  'X4. = X2*X7 + X3*X6 + k54*X5 - X4*X6 - X4*X7 - k45*X4',
               ##  'X5. = k45*X4 - k54*X5 - b5*X5'
               ##  ],
               ## ['X1. = a1 + 2*k21*X2 - 2*k12*X1^2',
               ##  'X2. = k12*X1^2 + X3*X8 + X5*X7 - k21*X2 - X2*X7 - X2*X8',
               ##  'X3. = X2*X7 + X4*X8 - X3*X7 - X3*X8',
               ##  'X4. = X5*X8 + X3*X7 - X4*X7 - X4*X8 ',
               ##  'X5. = X2*X8 + X4*X7 + k65*X6 - X5*X7 - X5*X8 - k56*X5',
               ##  'X6. = k56*X5 - k65*X6 - b6*X6'
               ##  ],
               ## ['X1. = a1 - k12*X1',
               ##  'X2. = k12*X1 + X3*X8 + X5*X7 - X2*X7 - X2*X8',
               ##  'X3. = X2*X7 + X4*X8 - X3*X7 - X3*X8',
               ##  'X4. = X5*X8 + X3*X7 - X4*X7 - X4*X8 ',
               ##  'X5. = X2*X8 + X4*X7 - X5*X7 - X5*X8 - k56*X5',
               ##  'X6. = k56*X5 - b6*X6'
               ##  ]
               ]



for f in system_sets:
    eq = dspace.Equations(f)
    # The next step is to construct an Equations object using this list of strings.
    ## Construct the System Design Space - with and without resolving cycles.
    # The previous system of equations can be analyzed using system design space.
    # First, we will define a system design space in the way we have defined it
    # in the previous examples, without resolving any cycles.
    
    ds_nocycles = dspace.DesignSpace(eq, name='Without Cycles')
    
    # Then, we will resolve the cycles in a second design space object and we will 
    # compare it with the design space where the cycles are not resolved. To resolve
    # cycles, the design space is constructed with the 'resolve_cycles=true' 
    # keyword argument.
    
    ds_cycles = dspace.DesignSpace(eq, resolve_cycles=True, resolve_codominance=True, name='With Cycles')#,constraints=constraints)
    
    ## Test for Validity
    # We calculate the valid cases without and with resolving cycles. As we can see
    # the valid cases are not the same for the two methods.  The set of valid cases
    # for the design space without resolving cycles is a subset of those found when
    # resolving cycles. 
    
    valid_no_cycles = ds_nocycles.valid_cases()
    valid_cycles = ds_cycles.valid_cases()
    
    for i in valid_cycles:
        c = ds_cycles(i)
        if '_' in i:
            break
    pvals=c.valid_interior_parameter_set()

    ds_cycles.draw_2D_slice_interactive(pvals,
                                        'X7',
                                        'X8',
                                        [1e-3, 1e3],
                                        [1e-3, 1e3],
                                        {i:[1e-5, 1e5] for i in ds_cycles.independent_variables if i not in ['X6', 'X7']},
                                        )
    
    fig = figure(1)
    clf()
    
    # We plot the design space with cycles first to get the colors of all the
    # cases and subcases. These colors are used to make the plot without resolving
    # the cycles to match the colorbars.
    
    ax = fig.add_subplot('122')       # We define the plot axes.
    colors = ds_cycles.draw_2D_slice(ax,
                             pvals,             # Pass the reference parameter set.
                             'X7',              # First, indicate the x-axis variable.
                             'X8',              # Second, indicate the y-axis variable.
                             [1e-3, 1e3],       # The range on the x-axis.
                             [1e-3, 1e3],       # The range on the y-axis.
                             expand_cycles=False
                             )
    ax.set_title(ds_cycles.name)

    ax = fig.add_subplot('121')       # We define the plot axes.
    ds_nocycles.draw_2D_slice(ax,
                      pvals,             # Pass the reference parameter set.
                      'X7',              # First, indicate the x-axis variable.
                      'X8',              # Second, indicate the y-axis variable.
                      [1e-3, 1e3],       # The range on the x-axis.
                      [1e-3, 1e3],       # The range on the y-axis.
                      color_dict=colors,
                      )
    ax.set_title(ds_nocycles.name)
    
    show()