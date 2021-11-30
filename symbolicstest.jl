import Pkg
Pkg.activate(".")
Pkg.add(["DifferentialEquations", "Plots",
        "ModelingToolkit","Latexify"])

using DifferentialEquations
using Plots
using Markdown
using ModelingToolkit
using Latexify

@parameters t α β γ δ
@variables x(t) y(t)
D = Differential(t)

eqs = [D(x) ~ α*x - β*x*y,
       D(y) ~ -δ*y + γ*x*y]

sys = ODESystem(eqs)

u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,0.75,3.0,0.75]

prob = ODEProblem(sys,u0,tspan,p)
sol = solve(prob)

plot(sol)