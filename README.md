This solves the Klein-Gordon equation $\frac{\partial^2 \phi}{\partial t^2} = \nabla \phi + V(\phi)$, based on the wave example in the [AMReX guided tutorials](https://github.com/AMReX-Codes/amrex-tutorials). 

The potential is $m^2 \phi^2$, and in the case of multiple scalar fields this is: $\phi^2 = \phi_1^2 + \phi_2^2 + ...$

The initial condition is a Gaussian pulse: $1.0 + A * exp(-r^2/\sigma)$


There are several properties of the scalar fields that are set in the input parameter file:
* wave.initial_amplitude - the initial amplitude of the Gaussian pulse, there should nfield integers here, separated by a space
* wave.initial_width - intial width of Gaussian pulse, again one value per scalar field
* wave.scalar_mass - scalar field mass, one per scalar field
* wave.tagging_criterion - value (in first derivative) at which cells will be tagged. 
* wave.v - verbosity level. Setting > 1 will output Conduit grid information and Conduit Blueprint files. 

For in-situ visualizations, you will also need to set:
* catalyst_input_script - path to the Python script with ParaView commands
* catalyst_options (optional) - pass optional flags to ParaView e.g. --enable-live for ParaView Live

The number of fields simulated is set in a macro NFIELDS defined in AmrLevelWave.H

Build Instructions: 
* Clone the AMReX repository [here](https://github.com/AMReX-Codes/amrex). AMReX has the following requirements: GNU make >=3.81, Python >=2.7, C++ compiler with C++17 support, Fortran compiler with Fortran 2003 standard support
* Set `AMREX_HOME` to the location of the AMReX directory
* Clone this repository. For the in-situ component, you will also need ParaView >= 5.9, Conduit, Catalyst 2.0 and make sure that `PARAVIEW_DIR` and `CONDUIT_DIR` are set to their install directories. 
* In the directory for this repo, edit the GNUmakefile for your system - here you can choose your compiler by changing the COMP variable (either gnu, intel or intel-llvm). You can also also change `USE_CATALYST=TRUE` for in-situ visualizations. If you want to debug by outputting Conduit Blueprint files you will also need to set `USE_CONDUIT= TRUE`. (Also set the verbosity level.)
* If `USE_CATALYST=TRUE`, you will need to define `$PARAVIEW_DIR` such athat it points to the directory where ParaView and Catalyst are installed, i.e. `$PARAVIEW_DIR/include/Catalyst-2.0` should contain the Catalyst2 header files. 
* Then run `make` 


From the original README: 

>The Laplacian operator is discretized dimension by dimension with a fourth-order stencil,

>$$\frac{\partial^2 u}{\partial x^2} = \left(-\frac{5}{2} u_i + \frac{4}{3} (u_{i-1} + u_{i+1}) - \frac{1}{12} (u_{i-2} + u_{i+2})\right) / \Delta x^2$$

>The time stepping is done with a Runge-Kutta method (RK2, RK3 or RK4).  
>In this test, the displacement at the x-direction boundaries is zero, and the it's periodic in the y-direction.  
>Note that refluxing is not implemented in this test code.
