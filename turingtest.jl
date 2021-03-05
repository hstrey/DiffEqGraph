using Turing
using Distributions
using Random
using Statistics

@model g(rn,T) = begin
    sigma ~ Uniform(0.0,20.0)
    
    rn ~ MvNormal(zeros(T),sigma)
end

N = 10
d = Normal(0.0, i)
data=rand(d,N)
chn = sample(g(data,N),NUTS(0.65),1000)

println("Std of chain: ",mean(Array(chn[:sigma])),"Std of data: ",std(data))
#println(macroexpand(g))
