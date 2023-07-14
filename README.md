# ETD.jl
Exploratory analysis of exponential time differencing

This program is based on a joint COMS 7300 final project (Spring 2017) by Toheeb Biala and Ashlin Harris.
It uses the exponential time differencing (ETD) scheme to obtain 
solutions to the fractional reaction diffusion equation.
The space derivatives are discretized using the matrix Transfer Technique
developed by  Ilic et al.[^1]. In time, we use the ETD schemes given by Kleefeld, Khaliq, and Wade[^2].

[^1]: Ilic, M.; Liu, F.; Turner, I. & Anh, V. Numerical approximation of a fractional-in-space 
    diffusion equation (II)- with nonhomogenous boundary conditions Fractional Calculus and 
    Applied Analysis, 2006.
[^2]: Kleefeld, B.; Khaliq, A. & Wade, B. An ETD Crank-Nicolson method for reaction-diffusion 
    systems Numerical Methods for Partial Differential Equations, Wiley Online Library, 2012,
    28, 1309-1335

