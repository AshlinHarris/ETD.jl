# ETD.jl
Exploratory analysis of exponential time differencing

This program is based on a joint COMS 7300 final project (Spring 2017) by Toheeb Biala and Ashlin Harris.
It uses the exponential time differencing (ETD) scheme to obtain 
solutions to the fractional reaction diffusion equation.
Space derivatives are discretized using the Matrix Transfer Technique (MTT)
developed by  Ilic et al.[^1].
For the time dimension, we use the ETD schemes given by Kleefeld, Khaliq, and Wade[^2].

[^1]: Ilic, M.; Liu, F.; Turner, I. & Anh, V. "Numerical approximation of a fractional-in-space 
    diffusion equation (II) - with nonhomogenous boundary conditions" _Fractional Calculus and 
    Applied Analysis_, 2006.
[^2]: Kleefeld, B.; Khaliq, A. & Wade, B. "An ETD Crank-Nicolson method for reaction-diffusion 
    systems" _Numerical Methods for Partial Differential Equations_, Wiley Online Library, 2012,
    28, 1309-1335

