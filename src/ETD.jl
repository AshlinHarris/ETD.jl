@info "Loading modules..."
using LinearAlgebra
using Plots

function main()
	# Space dimension
	Δx = 1
	X = -100:Δx:100
	N = length(X)

	# Time dimension
	Δt = 0.1
	T = 0:Δt:30
	M = length(T)

	Alphas = [2.0 1.99 1.9 1.7]
	final_slices = zeros(N,length(Alphas))
	furthest_spread = zeros(M,length(Alphas))
	for i in 1:length(Alphas)
		α = Alphas[i]
		@info "α = $α"

		u = initial_conditions(N,M)
		A = transfer_matrix(α, N, Δx)

		ETDBE!(u, N, M, Δt, A)
		savefig(figure_1(X, T, u, "Backward Euler, α = $α"), "surfaceBE$α.png")
		savefig(figure_2(X, T, u, "Backward Euler, α = $α"), "contourBE$α.png")
		
		ETDCN!(u, N, M, Δt, A)
		savefig(figure_1(X, T, u, "Crank Nicolson, α = $α"), "surfaceCN$α.png")
		savefig(figure_2(X, T, u, "Crank Nicolson, α = $α"), "contour$α.png")

		final_slices[:,i] = u[:,end]
		furthest_spread[:,i] = 1:M .|> y -> findlast(x->x>1E-6, u[:,y])
	end

 	savefig(plot(
		X, final_slices,
		label=Alphas .|> x->"α = $x",
		plot_title="density at t=30", #todo: fix hard-coded value
		xguide="Space (x)",
		yguide="Density (u)",
		#linecolor=:gnuplot,
		linewidth=3.0,
		levels=0.1:0.1:1.0,
		#zticks=0.1:0.1:1.0,
		#grid=false,
		gridalpha=1.0,
		gridlinewidth=0.30,
		minorgrid=true,
		minorgridalpha=1.0,
		minorgridlinewidth=0.10,
	), "final_slices.png")

 	savefig(plot(
		T, furthest_spread,
		label=Alphas .|> x->"α = $x",
		plot_title="Furthest spread",
		xguide="Time (t)",
		yguide="Space (x)",
		#linecolor=:gnuplot,
		linewidth=3.0,
		levels=0.1:0.1:1.0,
		#zticks=0.1:0.1:1.0,
		#grid=false,
		gridalpha=1.0,
		gridlinewidth=0.30,
		minorgrid=true,
		minorgridalpha=1.0,
		minorgridlinewidth=0.10,
	), "furthest_spread.png")

	return
end

function initial_conditions(N,M)
	u = zeros(N,M)

	midpoint = ceil(Int, N/2)
	bias = N ÷ 16
	u[(midpoint-bias):(midpoint+bias), 1].=0.8

	return u
end

function transfer_matrix(α, N, Δx)
	D = 0.1 # diffusion coefficient
	length1 = 1:1.0:N-2
	A1 = diagm(-D/Δx.^α*(2sin.(length1.*π/(2*(N-1)))).^α)
	P = sin.(length1'.*length1.*π/(N-1));
	return P\A1*P
end

# Reaction term
function f(v,t)
	K = 1.0
	r = 0.25
	return r*v.*(1 .- v/K)
end

# Exponential time differencing, Backward Euler Scheme
function ETDBE!(u, N, M, Δt, A)
	for j in 1:M-1
		u[2:N-1,j+1] = (I - Δt*A) \ (u[2:N-1,j] + Δt * f(u[2:N-1,j],j))
	end
end

# Exponential time differencing, Crank-Nicolson scheme
function ETDCN!(u, N, M, Δt, A)
	b = zeros(N-2,M-1)
	Nb = zeros(N-2,M-1)
	Nu = zeros(N-2,M-1)
	for j in 1:M-1
		Nb[:,j] = (2I - Δt*A) \ (4u[2:N-1,j] + 2Δt*f(u[2:N-1,j], j))
		b[:,j] = Nb[:,j] - u[2:N-1,j]
		Nu[:,j] =  (2I - Δt*A)\(Δt*(f(b[:,j],j+1) - f(u[2:N-1,j], j)))
		u[2:N-1,j+1] = b[:,j] + Nu[:,j]
	end
end

# Surface plot
function figure_1(X,T,u, title)
	return plot(
		X[1:end÷41:end],
		T[1:end÷21:end],
		u'[1:end÷21:end,1:end÷41:end],
		plot_title=title,
		xguide="Space (x)",
		yguide="Time (t)",
		zguide="Density (u)",
		linetype=:wireframe,
		#grid=false,
		gridalpha=1.0,
		gridlinewidth=0.30,
		minorgrid=true,
		minorgridalpha=1.0,
		minorgridlinewidth=0.10,
	)
end

# Contour plot
function figure_2(X,T,u, title)
	return contour(
		X,
		T,
		u',
		plot_title=title,
		xguide="Space (x)",
		yguide="Time (t)",
		zguide="Density (u)",
		color=:gnuplot,
		linewidth=3.0,
		levels=0.1:0.1:1.0,
		#zticks=0.1:0.1:1.0,
		#grid=false,
		gridalpha=1.0,
		gridlinewidth=0.30,
		minorgrid=true,
		minorgridalpha=1.0,
		minorgridlinewidth=0.10,
	)
end

main()

