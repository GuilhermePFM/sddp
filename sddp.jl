using JuMP, Clp
include("types.jl") 

function subproblem(states::States, FCF, a::Float64, t::Int, phase::Phase)
    m = Model(solver = ClpSolver())

    # Variables
    @variable(m, 0 <= g[1:2] <= 100)
    @variable(m, 0 <= u <= 150)
    @variable(m, 0 <= v <= 200)
    @variable(m, alpha >= 0)
    @variable(m, deficit >= 0)

    # -- Constraints
    # load balance
    rho = 1
    @constraint(m, loadbalance, g[1] + g[2] + u * rho + deficit == 150)
    
    # water balance
    s = 0
    @constraint(m, waterbalance, v + u + s == states.v[t-1] + a)

    # FCF cuts
    add_cuts!(m, FCF, t)

    # Objective function
    @objective(m, Min, 100*g[1] + 1000*g[2] + 100000 * deficit + alpha)

    # solve model
    status = solve(m)
    println("=== subproblema ==")
    println("phase = $phase")
    println("stage = $t")
    println("G = $(getvalue(g))")
    println("u = $(getvalue(u))")
    println("v = $(getvalue(v))")

    if phase == BACKWARD
        # return cut parameters
        γ1 = getdual(waterbalance)
        γ0 = getobjectivevalue(m)
        
        return γ1, γ0 
    elseif phase == SIMULATION
        # return model
        return m
    elseif phase == FORWARD
        # update state
        states.v[t] = getvalue(v)
    end

    return getobjectivevalue(m)
end

function add_cuts!(m::JuMP.Model, FCF, t::Int)
    # get cuts for state t
    cuts = FCF[t+1]
    nCuts = length(cuts)

    # get variables from problem
    alpha = getvariable(m, :alpha)
    v = getvariable(m, :v)

    # add the cuts
    @constraintref FCF_cstr[1:nCuts]
    for (i, (γ1, γ0)) in enumerate(cuts)
        @constraint(m, FCF_cstr[i], alpha >= γ0 + γ1 * v)
    end    
end

function forward!(prb::Problem, states::States, FCF, a_matrix::Array{Array{Float64,1},1})
    nStages = prb.nStages
    nForwards = prb.nForwards
    C_vec = zeros(nForwards)

    for t in 2:nStages
        for f in 1:nForwards
            a = a_matrix[t][f]
            C_vec[f] += subproblem(states, FCF, a, t, FORWARD)
        end
    end

    # build test statistics 
    Csup = sum(1/nForwards * C_vec)
    σ2 = std(C_vec)

    return Csup, σ2
end

function simulation!(prb::Problem, states::States, FCF, a_matrix::Array{Array{Float64,1},1})
    nStages = prb.nStages
    nForwards = prb.nForwards
    m_vec = []
    for t in 2:nStages
        f = 1
        a = a_matrix[t][f]
        m = subproblem(states, FCF, a, t, SIMULATION)
        push!(m_vec, m)
    end

    # print(m_vec[1])
end

function backward!(prb::Problem, states::States, FCF, a_matrix::Array{Array{Float64,1},1})
    nStages = prb.nStages
    nBackwards = prb.nBackwards

    Cinf = 0
    for i in 0:(nStages - 2)
        t = nStages - i
        vec_γ1 = []
        vec_γ0 = []

        # calculates FCF approximations
        for b in 1:nBackwards
            a = a_matrix[t][1]
            γ1, γ0 = subproblem(states, FCF, a, t, BACKWARD)

          push!(vec_γ1, γ1)
          push!(vec_γ0, γ0)
        end

        # updates FCF with expected value for approximations
        avrg_γ1 = sum(1/nBackwards * vec_γ1)
        avrg_γ0 = sum(1/nBackwards * vec_γ0) - avrg_γ1 * states.v[t-1]
        push!(FCF[t], (avrg_γ1, avrg_γ0))

        # 2 => first stage + 1 => FCF relative to first stage
        if t == 2
            # update Cinf
            Cinf = avrg_γ0
        end
    end

    return Cinf
end

function sddp()
    tolerance = 0.05
    nStages = 3  # sempre um a mais
    nForwards = 3
    nBackwards = 3
    prb = Problem(nStages, nBackwards, nForwards)
    
    # states begin at stage 2, to avoid 0 index error
    states = States(prb.nStages)
    # initial storage
    states.v[1] = 50

    # a_matrix = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 2.0, 3.0], [3.0, 3.0, 3.0]]
    inf_forward = [[0.0, 0.0, 0.0], [90, 100, 110], [10, 0, 0]]
    inf_backward = [[0.0, 0.0, 0.0], [50, 100, 150], [100, 50, 0]]
    inflw = Inflows(inf_forward, inf_backward)
    
    # FCF considers the FCF from last stage + 1
    nFCF = (prb.nStages + 1)
    FCF = [Tuple{Float64,Float64}[] for i in 1:nFCF]
    
    max_it = 10
    for i in 1:max_it
        Csup, σ2 = forward!(prb, states, FCF, inflw.forward)
        Cinf = backward!(prb, states, FCF, inflw.backward)

        # build confidence 95% interval
        # IC_ub = Csup + Z * sqrt(σ2 / prb.nForwards)
        # IC_lb = Csup - Z * sqrt(σ2 / prb.nForwards)
        IC_ub = Csup + 2 * σ2
        IC_lb = Csup - 2 * σ2

        # test for convergence
        if Cinf <= IC_ub &&  IC_lb <= Cinf
            println("Solution found!")
            println("Zsup = $(Csup)")
            println("Zinf = $(Cinf)")
            break
        end
    end

    # simulate 
    simulation!(prb, states, FCF, inflw.forward)
end