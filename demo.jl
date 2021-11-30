### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b55e8366-6fcd-11eb-11fd-77cc5e771425
begin
	import Pkg
	Pkg.activate(".")
	Pkg.add(["DifferentialEquations", "Plots", "PlutoUI","GraphPlot","LightGraphs","GraphRecipes","Markdown"])

	using DifferentialEquations
	using PlutoUI
	using Plots
	using GraphPlot
	using LightGraphs
	using GraphRecipes
	using Markdown
end

# ╔═╡ a2d4cc0a-75ea-11eb-02dc-53ffb118bfec
md"""
# Solving a simple system of differential equations
Here we are looking at two interacting variables x and y
Their change in time is governed by the interactions between x and y
and depends on 4 parameters: α, β, δ, γ
"""

# ╔═╡ 89659ab0-3f21-11eb-0eda-e99c5647e6f6
function lotka_volterra(du,u,p,t)
	  x, y = u
	  α, β, δ, γ = p
	  du[1] = dx = α*x - β*x*y
	  du[2] = dy = -δ*y + γ*x*y
end

# ╔═╡ fa33bece-3f22-11eb-2670-9fa1678c1504
u0 = [1.0,1.0]

# ╔═╡ 13b2b436-3f23-11eb-1059-9310363dceeb
tspan = (0.0,10.0)

# ╔═╡ 103757ea-3f22-11eb-0d50-f1434cabf86d
@bind alpha html"alpha <input type=range min=1.0 max=2.0 step=1e-2>"

# ╔═╡ 26b0c852-3f23-11eb-2d57-1da8f1831517
@bind beta html"beta <input type=range min=0.5 max=1.0 step=1e-2>"

# ╔═╡ 52e96334-3f23-11eb-3b16-0bab6a008978
@bind delta html"delta <input type=range min=2.0 max=4.0 step=1e-2>"

# ╔═╡ 7e08dc2a-3f23-11eb-11f3-6fab607a18c9
@bind gamma html"gamma <input type=range min=0.5.0 max=1.5 step=1e-2>"

# ╔═╡ 13b2eed8-3f23-11eb-2401-bd8e645e10bb
p = [alpha,beta,delta,gamma]

# ╔═╡ 13b3466c-3f23-11eb-33b5-6b885f523bcb
prob = ODEProblem(lotka_volterra,u0,tspan,p)

# ╔═╡ f2a32786-3f21-11eb-1a6f-1d3939e83184
sol = solve(prob)

# ╔═╡ 9c69556e-3f23-11eb-1f68-51527dd554e9
md"alpha $alpha, beta $beta, delta $delta gamma $gamma"

# ╔═╡ 044dd4fe-3f22-11eb-2b11-8bbb13cc851b
plot(sol,size=(3000,2000),thickness_scaling = 4)

# ╔═╡ e4e37a7a-75e9-11eb-2cbd-fb18edac47d5
md"""
# Breaking up the pieces
this is a very simple demonstration of how NeuroBlox would work under the hood
Our system of differential equations has two variables/nodes: x and y
These nodes interact with each other through the parameters: α, β, δ, γ
Graphically this could be expressed as arrows
"""

# ╔═╡ 82411708-75ec-11eb-0ba3-4b96c34c66c6
begin

g = [1 1;
     1 1]

graphplot(DiGraph(g), self_edge_size=0.15,nodeshape=:circle,names=["x","y"],size=(3000,2000),nodesize=1,thickness_scaling = 10)
end

# ╔═╡ 71aec586-75ef-11eb-1e2a-2d1921b0b80b
md"""
we have to now define the differential equations in terms of how the nodes interact
This can be done by the order of interaction.
For example: dx = α x + β x y
The first term is linear, the second is quadratic.  We can express this as matrixes.
"""

# ╔═╡ 22affc88-75f5-11eb-2ba3-87e944044a13
linear_matrix = ["α" ""
				 "" "δ"]

# ╔═╡ dd5d03c0-7628-11eb-2dc6-276df5b005ae
quadratic_matrix = ["β","γ"]

# ╔═╡ 79199746-75f5-11eb-3afd-91d158807973
nodes = String["x","y"]

# ╔═╡ 8eb417b6-75f5-11eb-3775-9b2c8795b5aa
function make_linear(lin_mat::Array{String,2}, n::Array{String,1})
	expressions = String[]
	for row in eachrow(lin_mat)
		expression = ""
		@show row
		for (element,item) in zip(row,n)
			@show element
			@show item
			if element != ""
				if expression == ""
					expression = element*"*"*item
				else
					expression = expression*"+"*element*"*"*item
				end
			end
		end
	push!(expressions, expression)
	end
	println("********")
	return expressions
end

# ╔═╡ 05b6da98-7629-11eb-3e66-a9fc58f7b805
function make_quad(quad_mat::Array{String,1}, n::Array{String,1})
	expressions = String[]
	for element in quad_mat
		expression = element
		if element!=""
			@show element
			for item in n
				@show item
				expression = expression*"*"*item
			end
		end
	push!(expressions, expression)
	end
	println("********")
	return expressions
end

# ╔═╡ 23a5e782-75f6-11eb-3d5d-bde9b5f6545c
make_linear(linear_matrix,nodes)

# ╔═╡ 98d51aea-7629-11eb-0d91-89bc7097fca5
make_quad(quadratic_matrix,nodes)

# ╔═╡ b901be18-7629-11eb-0783-b905e681ab8f
function make_equations(
		lin_mat::Array{String,2},
		quad_mat::Array{String,1},
		n::Array{String,1})
	linear = make_linear(lin_mat,n)
	quad = make_quad(quad_mat,n)
	expressions = String[]
	for (l,q,v) in zip(linear,quad,n)
		expression = "d"*v*" = "*l*" + "*q
		push!(expressions,expression)
	end
	return expressions
end
	

# ╔═╡ 6e9b82a4-762a-11eb-3a54-c5223231e0e2
make_equations(linear_matrix,quadratic_matrix,nodes)

# ╔═╡ 2181b2c2-76c0-11eb-1c96-69bd51f28fc1
para = ["α","β","δ","γ"]

# ╔═╡ d02be14e-7659-11eb-26a2-554021f09069
function generate_codestring(lin_mat,quad_mat,n,para)
	codestring = "(du,u,p,t)->begin "
	equations = make_equations(lin_mat,quad_mat,n)
	for nn in n
		codestring = codestring*nn*","
	end
	codestring = chop(codestring,tail=1)*"=u;"
	for p in para
		codestring = codestring*p*","
	end
	codestring = chop(codestring,tail=1)*"=p;"
	for (i,e) in enumerate(equations)
		codestring = codestring*"du["*string(i)*"]="*e*";"
	end
	codestring = codestring*"end"
	return codestring
end

# ╔═╡ 38c49f4a-76bc-11eb-048b-1902429bbb32
code = generate_codestring(linear_matrix,quadratic_matrix,nodes,para)

# ╔═╡ 80489b4c-76be-11eb-123e-7f13602eadc5
@eval f = $(Meta.parse(code))

# ╔═╡ e3609522-76be-11eb-28ae-ad88f00a1a5f
p2 = [alpha,-beta,-delta,gamma]

# ╔═╡ 96e66456-76be-11eb-1061-dba11e64dc63
prob2 = ODEProblem(f,u0,tspan,p2)

# ╔═╡ b2f5641e-76be-11eb-0f44-b91ee0cfc6b7
sol2 = solve(prob2)

# ╔═╡ cb711a9a-76be-11eb-175c-15cc343bdeef
plot(sol2)

# ╔═╡ Cell order:
# ╠═b55e8366-6fcd-11eb-11fd-77cc5e771425
# ╟─a2d4cc0a-75ea-11eb-02dc-53ffb118bfec
# ╠═89659ab0-3f21-11eb-0eda-e99c5647e6f6
# ╠═fa33bece-3f22-11eb-2670-9fa1678c1504
# ╠═13b2b436-3f23-11eb-1059-9310363dceeb
# ╠═13b2eed8-3f23-11eb-2401-bd8e645e10bb
# ╠═13b3466c-3f23-11eb-33b5-6b885f523bcb
# ╠═f2a32786-3f21-11eb-1a6f-1d3939e83184
# ╟─103757ea-3f22-11eb-0d50-f1434cabf86d
# ╟─26b0c852-3f23-11eb-2d57-1da8f1831517
# ╟─52e96334-3f23-11eb-3b16-0bab6a008978
# ╟─7e08dc2a-3f23-11eb-11f3-6fab607a18c9
# ╟─9c69556e-3f23-11eb-1f68-51527dd554e9
# ╠═044dd4fe-3f22-11eb-2b11-8bbb13cc851b
# ╠═e4e37a7a-75e9-11eb-2cbd-fb18edac47d5
# ╠═82411708-75ec-11eb-0ba3-4b96c34c66c6
# ╠═71aec586-75ef-11eb-1e2a-2d1921b0b80b
# ╠═22affc88-75f5-11eb-2ba3-87e944044a13
# ╠═dd5d03c0-7628-11eb-2dc6-276df5b005ae
# ╠═79199746-75f5-11eb-3afd-91d158807973
# ╠═8eb417b6-75f5-11eb-3775-9b2c8795b5aa
# ╠═05b6da98-7629-11eb-3e66-a9fc58f7b805
# ╠═23a5e782-75f6-11eb-3d5d-bde9b5f6545c
# ╠═98d51aea-7629-11eb-0d91-89bc7097fca5
# ╠═b901be18-7629-11eb-0783-b905e681ab8f
# ╠═6e9b82a4-762a-11eb-3a54-c5223231e0e2
# ╠═2181b2c2-76c0-11eb-1c96-69bd51f28fc1
# ╠═d02be14e-7659-11eb-26a2-554021f09069
# ╠═38c49f4a-76bc-11eb-048b-1902429bbb32
# ╠═80489b4c-76be-11eb-123e-7f13602eadc5
# ╠═e3609522-76be-11eb-28ae-ad88f00a1a5f
# ╠═96e66456-76be-11eb-1061-dba11e64dc63
# ╠═b2f5641e-76be-11eb-0f44-b91ee0cfc6b7
# ╠═cb711a9a-76be-11eb-175c-15cc343bdeef
