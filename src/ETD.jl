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

@info "Loading modules..."
using LinearAlgebra
using Plots

function figure_1(X,T,u, title)
	return plot(
		X[1:end÷41:end],
		T[1:end÷21:end],
		u'[1:end÷21:end,1:end÷41:end],
		linetype=:wireframe,
		plot_title=title,
		#grid=false,
		gridalpha=1.0,
		gridlinewidth=0.30,
		minorgrid=true,
		minorgridalpha=1.0,
		minorgridlinewidth=0.10,
		xguide="Space (x)",
		yguide="Time (t)",
		zguide="Density (u)",
	)
end

function main()

	K = 1.0
	r = 0.25
	D = 0.1

	Δx = 1
	Δt = 0.1
	X = -100:Δx:100
	T = 0:Δt:30

	N = length(X)
	M = length(T)

	f(v,t) = r*v.*(1 .- v/K)

	for α in [2.0 1.8]
		@info "α = $α"

		# Transfer matrix
		length1 = 1:1.0:N-2
		A1 = diagm(-D/Δx.^α*(2sin.(length1.*π/(2*(N-1)))).^α)
		P = sin.(length1'.*length1.*π/(N-1));
		A = P\A1*P

		u = zeros(N,M)
		# Initial condition
		midpoint = ceil(Int, N/2)
		bias = N ÷ 16
		u[(midpoint-bias):(midpoint+bias) , 1].=0.8
		
		# ETD, Crank-Nicolson scheme
		b = zeros(N-2,M-1)
		Nb = zeros(N-2,M-1)
		Nu = zeros(N-2,M-1)
		for j in 1:M-1
			Nb[:,j] = (2I - Δt*A) \ (4u[2:N-1,j] + 2Δt*f(u[2:N-1,j], j))
			b[:,j] = Nb[:,j] - u[2:N-1,j]
			Nu[:,j] =  (2I - Δt*A)\(Δt*(f(b[:,j],j+1) - f(u[2:N-1,j], j)))
			u[2:N-1,j+1] = b[:,j] + Nu[:,j]
		end
		savefig(figure_1(X, T, u, "Crank Nicolson, α = $α"), "CN$α.png")

		# ETD, Backward Euler Scheme
		for j in 1:M-1
			u[2:N-1,j+1] = (I - Δt*A) \ (u[2:N-1,j] + Δt * f(u[2:N-1,j],j))
		end
		savefig(figure_1(X, T, u, "Backward Euler, α = $α"), "BE$α.png")

		#savefig(contour(
		#	T,
		#	X,
		#	u,
		#	camera=(-30,20),
		#	c=:Reds,
		#	), "BE$α.png")
		#savefig(plot(u[:,end]), "end_BE$α.png")

	end

end # main

main()
	
