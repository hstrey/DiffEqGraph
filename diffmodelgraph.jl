### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ bc5b20e0-7de5-11eb-1cf0-8d1c48956180
begin
	import Pkg
	Pkg.activate(".")
	Pkg.add(["DifferentialEquations", "Plots", 			          "PlutoUI","GraphPlot","LightGraphs","GraphRecipes",
			"Markdown","SimpleWeightedGraphs","ModelingToolkit"])

	using DifferentialEquations
	using PlutoUI
	using Plots
	using GraphPlot
	using LightGraphs
	using GraphRecipes
	using Markdown
	using SimpleWeightedGraphs
	using ModelingToolkit
end

# ╔═╡ e7057d90-7de5-11eb-3552-879f6b025983
begin
	graph = SimpleWeightedDiGraph(2)
	e1 = SimpleWeightedEdge(1, 2, 1)
	e2 = SimpleWeightedEdge(2, 1, 4)
	add_edge!(graph, e1)
	add_edge!(graph, e2)
	graphplot(graph,curvature_scalar=0.1,arrow=true)
end

# ╔═╡ 3830b28e-7de6-11eb-1851-8be0f8abb0d1
begin
# define the connectome graph

tng = SimpleWeightedDiGraph(4)

Dinhibs = [0.1, 0.1]
Dgaps = fill(-0.05, 1)

# Add mutual inhibition
for (ventral, dorsal, Dinhib) in zip(1:2, 3:4, Dinhibs)
    e1 = SimpleWeightedEdge(ventral, dorsal, Dinhib)
    e2 = SimpleWeightedEdge(dorsal, ventral, Dinhib)
    LightGraphs.add_edge!(tng, e1)
    LightGraphs.add_edge!(tng, e2)
end

# We handle gaps separately, since the head oscillator 
# does not have a forcing neuron pair upstream
for (ventral, dorsal, Dgap) in zip(2:2, 4:4, Dgaps)
    e1 = SimpleWeightedEdge(ventral - 1, ventral, Dgap)
    e2 = SimpleWeightedEdge(dorsal - 1, dorsal, Dgap)
    LightGraphs.add_edge!(tng, e1)
    LightGraphs.add_edge!(tng, e2)
end
	graphplot(tng, names=1:4, curvature_scalar=0.1,arrow=true)
end

# ╔═╡ ad8d8140-7ded-11eb-1ddd-abad58c0ec28
begin
# define the unit equation

#@parameters t g e b

@variables v(t) w(t) F(t)

#@derivatives D'~t

f(v) = v - v^3/3

fhn = [
    D(v) ~ f(v) - w + F,
    D(w) ~ e * (v - g * w + b)
]

end

# ╔═╡ bfc7c8ba-7df3-11eb-388b-5bc57d8b13b9
begin
function couple(system_to, systems_from, weights)
    return [0 ~ system_to.F +  sum(weights .* getproperty.(systems_from, :v))]
end

# generate a single equation spec from the graph

# first, create all relevant systems
systems = [
    ODESystem(
        fhn, t, [v,w,F], [g,e,b];
        name = Symbol(
                "n" * string(i) # ModelingToolkit.map_subscripts(string(i))
        )
    )
    for i in vertices(tng)
]
end

# ╔═╡ 9a69bc72-7df3-11eb-0bec-ff1b641a8b5a
begin
# populate couplings from edge

couplings = Equation[]

edgs = edges(tng)
graph_weights = sparse(weights(tng))
for vertex in vertices(tng)
    sys_to = systems[vertex]
    weights = []
    systems_from = ODESystem[]
    for neighbor in inneighbors(tng, vertex)
        push!(systems_from, systems[neighbor])
        push!(weights, graph_weights[neighbor, vertex])
    end
    append!(
        couplings,
        couple(
            systems[vertex],
            systems_from,
            weights
        )
    )
end

connected = ODESystem(
    couplings,
    t,
    [],
    [];
    systems = systems
)

# Construct the initial condition, special-casing the first
# 2 neurons which comprise the head oscillator

u0 = [
    systems[1].v => -1.0,
    systems[1].w => -0.51,
    systems[1].F => -0.01,
    systems[2].v =>  1.0,
    systems[2].w => -0.49,
    systems[2].F => 0.01,
]


for sys in systems[3:end]
    append!(u0,
        [
        sys.v =>  1.0,
        sys.w => -0.49,
        sys.F => 0.01,
        ]
    )
end

# Construct the parameters, which are the same across
# all systems in this case

p0 = Pair{Operation, Float64}[]

for sys in systems
    append!(
        p0,
        [
            sys.g => 0.8,
            sys.b => 0.46,
            sys.e => 0.04,
        ]
    )
end
end

# ╔═╡ bd69611a-7ded-11eb-2c74-2d05e6f0eb62
begin
	prob = ODEProblem(connected, u0, (0.0, 2000.0), p0)
	sol = solve(prob, Rodas5())

	plot(sol; vars = collect(1:3:12))
end

# ╔═╡ Cell order:
# ╠═bc5b20e0-7de5-11eb-1cf0-8d1c48956180
# ╠═e7057d90-7de5-11eb-3552-879f6b025983
# ╠═3830b28e-7de6-11eb-1851-8be0f8abb0d1
# ╠═ad8d8140-7ded-11eb-1ddd-abad58c0ec28
# ╠═bfc7c8ba-7df3-11eb-388b-5bc57d8b13b9
# ╠═9a69bc72-7df3-11eb-0bec-ff1b641a8b5a
# ╠═bd69611a-7ded-11eb-2c74-2d05e6f0eb62
