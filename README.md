Design Space Toolbox V2 (Python Package)
==============================================

Copyright (C) 2013-2015 Jason G. Lomnitz, Department of Biomedical Engineering,
University of California, Davis, California, United States. All rights reserved.
E-mail: <jlomn@ucdavis.edu>.

This is part of the Design Space Toolbox V2 project, a software implementation of the System Design Space method developed by Jason G. Lomnitz, originally developed in the laboratory of Michael A. Savageau (for examples of this methodology, see [1-6]). This method decomposes complex nonlinear systems into a finite number of tractable nonlinear subsystems.

The Design Space Toolbox V2 Python Interface is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Getting Started
---------------------

To use the Design Space Toolbox V2 from Python, you must have Python 2.7.x, as well as a working version of the **Design Space Toolbox V2 C Library** (an open-source library to be released early 2016). Access to the Design Space Toolbox V2 C Library available upon request by contacting Jason G. Lomnitz (jlomn@ucdavis.edu). In addition, plot generation and data visualization depends on the NumPy (www.numpy.org), Scipy (www.scipy.org) and Matplotlib (www.matplotlib.org) packages.

The preferred installation process is using the Design Space Toolbox V2 update and installation script (see https://bitbucket.org/jglomnitz/toolbox-update-script for detailed instructions). This script manages the Git repositories associated with the C Library, Python Interface, as well as a modified version of the GNU Linear Programming Kit (GLPK) library used by the Design Space Toolbox Project.

Examples of System Analysis Using the Design Space Toolbox
------------------------------------------------------------------------------------

The Design Space Toolbox V2 Python Package is a standard Python package that includes Python code and a Python-C bridge extension. The current version of the package lacks documentation; however, there are examples available that highlight the array of objects and methods for the analysis of non-linear systems.

These examples are found in the `dspace.examples` sub-package, and can be imported directly into an interactive python shell. There are a total of 5 examples that highlight the following aspects of the software package.

Examples:

1. Demonstrates the basic steps and functions to construct, analyze and explore a design space object using a very simple system.

2. Reiterates basic steps and demonstrates more advanced functions of the Design Space Toolbox.

3. Demonstrates the analysis of intersecting, or overlapping, cases, indicative of multiple steady states with potential for multiple dynamic behaviors.

4. Emphasizes the analysis of dominant sub-systems.

5. Demonstrates a higher-level api for analyzing systems in the design space framework.

To run an example, simply start the python interpreter and import the example.

E.g., for example 1:

    import dspace.examples.example_1

Once imported, you can obtain the path for the example's source by typing

    dspace.examples.example_1

at the command prompt.	


References
---------------

1. Savageau MA, Coelho PMBM, Fasani RA, Tolla DA, and Salvador A (2009) Phenotypes and tolerances in the design space of biochemical systems. _Proc. Natl. Acad. Sci. U.S.A._ 106, 6435–6440.

2. Savageau MA, and Lomnitz JG (2013) _Deconstructing Complex Nonlinear Models in System Design Space_, in Discrete and Topological Models in Molecular Biology (Jonoska, N., and Saito, M., Eds.). Springer.

4. Fasani RA, and Savageau MA (2010) Automated construction and analysis of the design space for biochemical systems. _Bioinformatics_ 26:2601–2609.

5. Lomnitz JG, and Savageau, MA (2013) Phenotypic deconstruction of gene circuitry. _Chaos_ 23, 025108.

6. Lomnitz JG, and Savageau MA (2014) Strategy Revealing Phenotypic Differences among Synthetic Oscillator Designs. _ACS Synth Biol_ 3(9):686–701.