using DifferentialEquations
using Plots

function lotka_volterra(du,u,p,t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α*x - β*x*y
    du[2] = dy = -δ*y + γ*x*y
end

u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,0.75,3.0,1.0]
prob = ODEProblem(lotka_volterra,u0,tspan,p)
sol = solve(prob)

s = plot(sol)

