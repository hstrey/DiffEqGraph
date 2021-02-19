### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b10b8954-3f20-11eb-0b6b-69b77f96f60c
using DifferentialEquations

# ╔═╡ 810e594c-3f21-11eb-331b-ef943eab7fba
using Plots

# ╔═╡ b55e8366-6fcd-11eb-11fd-77cc5e771425
import Pkg

# ╔═╡ abc02e66-6fd4-11eb-1c30-d121ede93969
Pkg.add("DifferentialEquations")

# ╔═╡ c89acb38-6fcd-11eb-3c3a-41c496b74d7e
Pkg.add("Plots")

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
plot(sol)

# ╔═╡ Cell order:
# ╠═b55e8366-6fcd-11eb-11fd-77cc5e771425
# ╠═abc02e66-6fd4-11eb-1c30-d121ede93969
# ╠═b10b8954-3f20-11eb-0b6b-69b77f96f60c
# ╠═c89acb38-6fcd-11eb-3c3a-41c496b74d7e
# ╠═810e594c-3f21-11eb-331b-ef943eab7fba
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
