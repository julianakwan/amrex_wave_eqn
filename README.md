This solves the Klein-Gordon equation $\frac{\partial^2 phi}{\partial t^2} = \nabla \phi$ + V(\phi), based on AMReX's wave example (see https://github.com/AMReX-Codes/amrex-tutorials). 

The initial condition are a Gaussian pulse with the amplitude and width as free parameters in the input parameter file. 

To run multiple scalar fields, use the nfields parameter in the input parameter file. You will also need to specify the initial width and amplitude for each. 

The potential is assumed to be $m^2 \phi^2$, and in the case of multiple scalar fields this is: $\phi^2 = \phi_1^2 + \phi_2^2 + ...$

m is also a free parameter in the input params file. 

This is from the original README: 
The Laplacian operator is discretized dimension by dimension with a fourth-order stencil,

$$\frac{\partial^2 u}{\partial x^2} = \left(-\frac{5}{2} u_i + \frac{4}{3} (u_{i-1} + u_{i+1}) - \frac{1}{12} (u_{i-2} + u_{i+2})\right) / \Delta x^2$$

The time stepping is done with a Runge-Kutta method (RK2, RK3 or RK4).  In
this test, the displacement at the x-direction boundaries is zero, and the
it's periodic in the y-direction.  Note that refluxing is not implemented in
this test code.
