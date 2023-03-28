using DifferentialEquations, LinearAlgebra
using Plots, LaTeXStrings


USE_TIKZ = true
if (USE_TIKZ)
    pgfplotsx()
    plt_ext = "tikz"
    plt_size= (350,250)
else
    #plotly()
    gr()
    plt_ext = "png"
    plt_size= (800,450)
end

include("aux.jl")



## First simulation

my_params = (   
    SIZE_q = 2, # number of populations - must be 2
    SIZE_n = (2,2), # number of strategies per population - a tuple of integers >0    
    θ = Float64(0.0002),
    γ = Float64(0.14),
    δ = Float64(0.0005),
    M = Float64.([1.0 0.3; 0.1 1.0]),
    #M = Float64.([1.0 0.2; 0.0 1.0]),
    β = [ Float64.([0.15;0.19]),Float64.([0.15;0.19]) ],
    r = [ Float64.([0.15;0.15]),Float64.([0.15]) ],
    c = [ Float64.([0.35;0.0]),Float64.([0.40;0.0]) ],
    υ = 4.0
    )


@assert all( all(diff(my_params.β[i]).>0)  for i=1:2) #betas elements are increasing
@assert all( all(diff(my_params.c[i]) .<0) for i=1:2)     #costs elements are increasing
@assert minimum(diag(diagm([my_params.β[i][1] for i=1:2])*my_params.M))>my_params.δ+my_params.θ+my_params.γ #assumption 1

#assumption 2
for i=1:2
    h(β,c,k) = (c[k]-c[k+1])/(β[k+1]-β[k])
    β = my_params.β[i]
    c = my_params.c[i]
    for j=1:(my_params.SIZE_n[i]-2)
        @assert h(β,c,j)>=h(β,c,j+1) "error on i=$(i), c: j=$(j)"
    end
end

I0   = Vector{Float64}([0.001; 0.0011])
x0_1 = Vector{Float64}([1.0;0.0])
x0_2 = Vector{Float64}([1.0;0.0])
q0   = Vector{Float64}([0.0;0.0])

cstar = 0.30


u0 = [I0,[x0_1,x0_2],q0];

my_params = (;my_params...,xstar=find_xstar(my_params,cstar))
my_params = (;my_params...,r=find_rΔ(my_params.xstar))
[ my_params.c[i]'my_params.xstar[i] for i=1:2 ]
#plot!()

dx0 = formatted2vec(u0)
df(dx0,formatted2vec(u0),my_params,0.0)

prob = ODEProblem(df,formatted2vec(u0),[0.0,1500],my_params)

#sol = solve(prob,Euler(), dt=0.1)
sol = solve(prob,AutoTsit5(Rosenbrock23()), save_everystep=true, saveat=1.0)

plot(
    plot(sol, idxs=[1,2], label=[L"I_1" L"I_2"], title="Infected", ylims=[0,0.40]),
    #plot(sol, idxs=[3,4,5,6], label=[L"x_1^1" L"x_2^1" L"x_1^2" L"x_2^2"], title="Social State, υ=$(my_params.υ)", legend_position=:outerright),
    plot(sol, idxs=[3,5], label=[L"x_1^1" L"x_1^2"], title="Social State, υ=$(my_params.υ)"),
    #plot(sol, idxs=[7,8], label=[L"q_1" L"q_2"], title=L"dynamic payoff state ($q$)"),
    size=plt_size,
    #layout = @layout [a ; b c]
    layout = @layout [a ; b]   
)

savefig("fig1_v3.$(plt_ext)")

ylims!(-1.0,1.0)
xlims!(2800.0,3000.0)











## Second simulation

my_params = (   
    SIZE_q = 2, # number of populations - must be 2
    SIZE_n = (2,1), # number of strategies per population - a tuple of integers >0    
    θ = Float64(0.0002),
    γ = Float64(0.14),
    δ = Float64(0.0005),
    M = Float64.([1.0 0.3; 0.1 1.0]),
    β = [ Float64.([0.15;0.19]),Float64.([0.19]) ],
    r = [ Float64.([0.15;0.15]),Float64.([0.0]) ],
    c = [ Float64.([0.35;0.0]),Float64.([0.0]) ],
    υ = 4.0
    )



I0   = Vector{Float64}([0.001;0.0011])
x0_1 = Vector{Float64}([1.0;0.0])
x0_2 = Vector{Float64}([1.0])
q0   = Vector{Float64}([0.0;0.0])

cstar = 0.30

u0 = [I0,[x0_1,x0_2],q0];

my_params = (;my_params...,xstar=find_xstar(my_params,cstar))


prob = ODEProblem(df,formatted2vec(u0),[0.0,2000.0],my_params)
sol = solve(prob,AutoTsit5(Rosenbrock23()), save_everystep=true, saveat=1.0)

plot(
    plot(sol, idxs=[1,2], label=[L"I_1" L"I_2"], title="Infected", ylims=[0,0.40]),
    #plot(sol, idxs=[3,4,5,6], label=[L"x_1^1" L"x_2^1" L"x_1^2" L"x_2^2"], title="Social State, υ=$(my_params.υ)", legend_position=:outerright),
    #plot(sol, idxs=[3,5], label=[L"x_1^1" L"x_1^2"], title="Social State, υ=$(my_params.υ)", legend_position=:outerright),
    plot(sol, idxs=[3], label=L"x_1^1", title="Social State, υ=$(my_params.υ)"),
    #plot(sol, idxs=[6], label=L"q_1", title=L"dynamic payoff state ($q$)"),
    size=plt_size,
    #layout = @layout [a ; b c]
    layout = @layout [a ; b]   
)

savefig("fig2_v3.$(plt_ext)")
