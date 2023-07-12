# Program based on joint COMS 7100 Final Project by Toheeb Biala and Ashlin Harris

# This program uses the exponential time differencing (ETD) scheme to obtain 
# solutions to the fractional reaction diffusion equation.
# The space derivatives are discretized using the matrix Transfer Technique
# developed by  [1]. In time, we use the ETD schemes given in [2].
# [1] Ilic, M.; Liu, F.; Turner, I. & Anh, V. Numerical approximation of a fractional-in-space 
#     diffusion equation (II)- with nonhomogenous boundary conditions Fractional Calculus and 
#     Applied Analysis, 2006.
# [2] Kleefeld, B.; Khaliq, A. & Wade, B. An ETD Crank-Nicolson method for reaction-diffusion 
#     systems Numerical Methods for Partial Differential Equations, Wiley Online Library, 2012,
#     28, 1309-1335

using LinearAlgebra
using Plots

K = 1
r = 0.25
D = 0.5

Δx = 0.1
Δt = 0.01
X = -50:Δx:50
T = 0:Δt:1

N = length(X)
M = length(T)

f(v,t) = r*v.*(1 .- v/K)

for α in [2.0 1.8]

	b = zeros(N-2,M-1)
	Nb = zeros(N-2,M-1)
	Nu = zeros(N-2,M-1)
	length1 = 1:1.0:N-2

	A1 = diagm(-D/Δx.^α*(2sin.(length1.*π/(2*(N-1)))).^α)
	P = sin.(length1'.*length1.*π/(N-1));
	#savefig(plot(P), "P$α.png")
	A = P\A1*P

	u = zeros(N,M)
	# Initial condition
	midpoint = ceil(Int, N/2)
	bias = N ÷ 16
	u[(midpoint-bias):(midpoint+bias) , 1].=0.8
	
	# ETDCN scheme
	for j=1:M-1
		Nb[:,j] = (2I - Δt*A) \ (4u[2:N-1,j] + 2Δt*f(u[2:N-1,j], j))
		b[:,j] = Nb[:,j] - u[2:N-1,j]
		Nu[:,j] =  (2I - Δt*A)\(Δt*(f(b[:,j],j+1) - f(u[2:N-1,j], j)))
		u[2:N-1,j+1] = b[:,j] + Nu[:,j]
	end
	savefig(surface(u, camera=(80,20)), "CN$α.png")
	savefig(plot(u[:,end]), "end_CN$α.png")

	# Backward Euler
	for j=1:M-1
		u[2:N-1,j+1] = (I - Δt*A)\(u[2:N-1,j] + Δt*f(u[2:N-1,j],j)); 

	end
	savefig(surface(u, camera=(80,20)), "BE$α.png")
	savefig(plot(u[:,end]), "end_BE$α.png")

end
	
