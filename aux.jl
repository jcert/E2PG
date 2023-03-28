using DifferentialEquations, LinearAlgebra
using JuMP, GLPK, Optim, NLsolve
using ForwardDiff
using Plots, LaTeXStrings

formatted2vec(u) = reduce(vcat,reduce(vcat,u))

function vec2formatted(p, u)


    I  = u[1:p.SIZE_q]
    x_1 = u[p.SIZE_q.+(1:p.SIZE_n[1])]
    x_2  = u[(p.SIZE_q+p.SIZE_n[1]).+(1:p.SIZE_n[2])]
    q  = u[(p.SIZE_q+sum(p.SIZE_n[1:2])).+(1:p.SIZE_q)]

    [I,[x_1,x_2],q]
end


function F(params,u,t,i,j)
    β = params.β
    r = params.r

    q  = u[3]    # dynamic payoff state

    p = [[0];[0]]
    p = [ β[1].*q[1]+r[1], β[2].*q[2]+r[2] ]

    return p[j][i]
end

function smith(c,u,t,i,k; λ=1.0, τ=1.0 )
    x =  u[2]
    
    mysum = 0.0
    for j ∈ 1:c.SIZE_n[k]
        mysum +=  x[k][j]*max(min(λ*(F(c,u,t,i,k)-F(c,u,t,j,k)), τ),0.0)
    end

    for j ∈ 1:c.SIZE_n[k]
        mysum += -x[k][i]*max(min(λ*(F(c,u,t,j,k)-F(c,u,t,i,k)), τ),0.0)
    end
    
    mysum
end

function find_xstar(params,cstar)

    B(x) = diagm([ params.β[i]'*x[i,1:params.SIZE_n[i]]  for i=1:params.SIZE_q])*params.M 

    function f(cstar,cstar_1)
        m = Model(GLPK.Optimizer)

        @variable(m, x[1:params.SIZE_q,1:maximum(params.SIZE_n)]>=0)
        for i=1:params.SIZE_q
            @constraint(m, sum(x[i,(params.SIZE_n[i]+1):end])==0)
        end
    
        for i=1:params.SIZE_q
            @constraint(m, sum(x[i,:])==1 ) 
        end
        
        @constraint(m, params.c[1]'*x[1,1:params.SIZE_n[1]]<=cstar_1)
        @constraint(m, params.c[2]'*x[2,1:params.SIZE_n[2]]<=cstar-cstar_1)
        
        @objective(m, Min, sum( params.β[i]'*x[i,1:params.SIZE_n[i]]  for i=1:params.SIZE_q))
    
        optimize!(m)

        value.(x)
    end


    T = 0.0:0.001:1.0
    aux = [ (maximum(abs.(eigvals(B(f(cstar,cstar_1))))),cstar_1) for cstar_1=T.*cstar] 
    best_cstar_1 = (T*cstar)[findmin(x->x[1], aux)[2]]
    plot(aux)
    [ f(cstar,best_cstar_1)[i,1:params.SIZE_n[i]] for i=1:params.SIZE_q]
end


function find_rΔ(xstar;ρ=1.0)
    
    rΔ = [[],[]]
    for i=1:2
        rΔ[i] = ρ.*(xstar[i].>0)
    end
    rΔ
    
end

function find_equilibrium_I(params,u)
    θ = params.θ
    γ = params.γ
    δ = params.δ
    M = params.M
    β = params.β
    
    B = diagm([ β[i]'*u[2][i]  for i=1:params.SIZE_q])*M 

    I = ones(params.SIZE_q)
   
    for k=1:50
        Q = (B-δ*LinearAlgebra.I)*I
        D = (B-δ*LinearAlgebra.I)*I+(θ+γ).*ones(params.SIZE_q)
    
        I = Q./D
         
        #@show I
    end
    
    #epidemics equilibrium
    I

end

function SIS_storage(params,u,xstar)
    θ = params.θ
    γ = params.γ
    δ = params.δ
    M = params.M
    β = params.β
    r = params.r
    υ = params.υ

    B(x) = diagm([ β[i]'*x[2][i]  for i=1:params.SIZE_q])*M 

    I  = u[1]    # epidemic state

    Bcurr= B(u)
    Itarget = find_equilibrium_I(params,u)

    w1 = Bcurr[2,1]
    w2 = Bcurr[1,2]

    w1*(I[1]-Itarget[1]-Itarget[1]*log(I[1]/Itarget[1]))
    +w2*(I[2]-Itarget[2]-Itarget[2]*log(I[2]/Itarget[2]))
    +υ*norm(B(u)-B(xstar))

end

#explicit
function ∂w_∂I_∂B_by_∂fi_ex(params,u,i)
    
    θ = params.θ
    γ = params.γ
    δ = params.δ
    M = params.M
    β = params.β
    r = params.r

    B(x) = diagm([ β[i]'*x[2][i]  for i=1:params.SIZE_q])*M 

    Bcurr = B(u) 

    Ieq = find_equilibrium_I(params, u)
    Seq = 1 .- Ieq


    ∂I = -diagm(Seq)*diagm(LinearAlgebra.I(params.SIZE_q)[i,:])*M*Ieq
    ∂I = ( -(γ+θ)*LinearAlgebra.I+diagm(Seq)*Bcurr-δ*diagm(Seq)-diagm((Bcurr-δ*LinearAlgebra.I)*Ieq)) \ ∂I 

    ∂w1 = -∂I[2]*Ieq[1]*Bcurr[2,1]+Seq[2]*∂I[1]*Bcurr[2,1]
    if i==1
        ∂w1 = ∂w1 + Ieq[1]*Seq[2]*M[2,1]
    end
        
    ∂w2 = -∂I[1]*Ieq[2]*Bcurr[1,2]+Seq[1]*∂I[2]*Bcurr[1,2]
    if i==2
        ∂w2 = ∂w2 + Ieq[2]*Seq[1]*M[1,2]
    end

    (∂w=[∂w1,∂w2],∂I=∂I)
end


B(params, x) = diagm([ params.β[i]'*x[2][i]  for i=1:params.SIZE_q])*params.M

function df(dv,v,params,t)
    u = vec2formatted(params,v)

    du = deepcopy(u)

    θ = params.θ
    γ = params.γ
    δ = params.δ
    M = params.M
    β = params.β
    υ = params.υ
    r = params.r
    xstar = params.xstar

    B(x) = diagm([ β[i]'*x[2][i]  for i=1:params.SIZE_q])*M 

    I  = u[1]    # epidemic state
    x1 = u[2][1] # social state - pop 1
    x2 = u[2][2] # social state - pop 2
    q  = u[3]    # dynamic payoff state


    #enforce set invariance
    I .= min.(max.(I,0.0),1.0) # guarantees that each component lies on [0,1]
    for i=1:params.SIZE_q
        u[2][i] .= min.(max.(u[2][i],0.0),1.0)  # guarantees that each component lies on [0,1]
        u[2][i] .= u[2][i]./sum(u[2][i])        # guarantees that the strategies of the population add up to 1
    end



    #epidemics
    du[1] = -(θ+γ).*I+diagm(1 .-I)*(B(u)*I-δ.*I)
    
    #evolutionary dynamics
    for i=1:params.SIZE_n[1]
        du[2][1][i] = smith(params,u,t,i,1; λ=0.1,τ=0.1 )
    end
    for i=1:params.SIZE_n[2]
        du[2][2][i] = smith(params,u,t,i,2; λ=0.1,τ=0.1 )
    end
        

    #dynamic payoff
    w = [ B(u)[2,1], B(u)[1,2]]
    ∂  = [ ∂w_∂I_∂B_by_∂fi_ex(params,u,1), ∂w_∂I_∂B_by_∂fi_ex(params,u,2) ]
    Ieq = find_equilibrium_I(params,u)
    #@show I, Ieq
    for i=1:params.SIZE_q
        if length(params.xstar[i])>1
            du[3][i] = -sum( ∂[i].∂w[k]*(I[k]-Ieq[k]-Ieq[k]*(log(I[k]/Ieq[k])))+(w[k]*log(Ieq[k]/I[k]))*∂[i].∂I[k] for k=1:2)
            du[3][i] = du[3][i] - 2*υ*β[i]'*(u[2][i]-xstar[i])*(M*M')[i,i]
        else
            du[3][i] = 0.0
        end
    end
    
    dv .= formatted2vec(du)
    nothing
end


