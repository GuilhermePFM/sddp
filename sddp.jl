using JuMP, Clp
include("types.jl") 
include("read_inputs.jl")

function subproblem(prb::Problem, states::States, FCF, a::Array{Float64,1}, dem, t::Int, f::Int, phase::Phase)
    gmax = prb.gmax
    gcost = prb.gcost
    umax = prb.umax
    vmax = prb.vmax
    
    m = Model(solver = ClpSolver())

    # Variables
    @variable(m, 0 <= g[i=1:prb.nGen] <= gmax[i])
    @variable(m, 0 <= u[i=1:prb.nHyd] <= umax[i])
    @variable(m, 0 <= v[i=1:prb.nHyd] <= vmax[i])
    @variable(m, alpha[i=1:prb.nHyd] >= 0)
    @variable(m, deficit >= 0)

    # -- Constraints
    # load balance
    rho = 1
    @constraint(m, loadbalance, sum(g) + sum(u) * rho + deficit == dem)
    
    # water balance
    s = 0
    @constraint(m, waterbalance[i=1:prb.nHyd], v[i] + u[i] + s == states.v[t-1,f,i] + a[i])

    # FCF cuts
    if t+1 <= length(states.v[:,1,1])
        add_cuts!(m, FCF, t)
    end

    # Objective function
    @objective(m, Min, sum(gcost.*g) + 100000 * deficit + sum(alpha))

    # solve model
    status = solve(m)
    
    if phase == BACKWARD
        # return cut parameters
        γ1 = getdual(waterbalance)
        γ0 = getobjectivevalue(m)
        
        return γ1, γ0 
    elseif phase == SIMULATION
        println("=== subproblema ==")
        println("phase = $phase")
        println("stage = $t")
        println("G = $(getvalue(g))")
        println("u = $(getvalue(u))")
        println("v = $(getvalue(v))")
        
        # return model
        return m
    elseif phase == FORWARD
        # update state
        states.v[t,f,:] = getvalue(v)
    end
    
    ICF = sum(gcost.*getvalue(g)) + 100000 * getvalue(deficit)
    return ICF
end

function add_cuts!(m::JuMP.Model, FCF, t::Int)
    # get cuts for state t
    cuts = FCF[t+1,1]
    nHyd = 1#length(cuts[1,:])
    nCuts = length(cuts)

    # get variables from problem
    alpha = getvariable(m, :alpha)
    v = getvariable(m, :v)

    # add the cuts
    @constraintref FCF_cstr[1:nCuts,1:nHyd]
    #    for h in 1:nHyd
    h=1

    for (i, (γ1, γ0)) in enumerate(cuts)
            @constraint(m, FCF_cstr[i,h], alpha[h] >= γ0 + γ1 * v[h])
    end    
    #       end

end

function forward!(prb::Problem, states::States, FCF, a_matrix::Array{Array{Float64,2},1})
    nStages = prb.nStages
    nForwards = prb.nForwards
    C_vec = zeros(nStages, nForwards)

    for t in 2:nStages
        for f in 1:nForwards
            a = a_matrix[t][f,:]
            dem = prb.dem[t][1]
            C_vec[t,f] += subproblem(prb, states, FCF, a, dem, t, f, FORWARD)
        end
    end
    # build test statistics 
    Csup = sum(1/nForwards * C_vec)
    Cinf = sum(1/nForwards * C_vec[2,:])
    σ2 = std(C_vec)

    return Csup, σ2, Cinf
end

function simulation!(prb::Problem, states::States, FCF, a_matrix::Array{Array{Float64,2},1})
    nStages = prb.nStages
    nForwards = prb.nForwards
    m_vec = []
    g1_vec = zeros(nStages)
    g2_vec = zeros(nStages)
    u_vec = zeros(nStages)
    v_vec = zeros(nStages)
    cparcial = zeros(nStages)
    
    for t in 2:nStages
        f = 1
        a = a_matrix[t][f,:]
        dem = prb.dem[t][1]
        m = subproblem(prb,states, FCF, a, dem, t, f, SIMULATION)
        if t==2
            writeLP(m, "etapa1.lp", genericnames=false)
        end
        if t==nStages
            writeLP(m, "etapa2.lp", genericnames=false)
        end
        push!(m_vec, m)
        g1_vec[t] = getvalue(getvariable(m, :g))[1]
        g2_vec[t] = getvalue(getvariable(m, :g))[2]
        u_vec[t] = getvalue(getvariable(m, :u))[1]
        v_vec[t] = getvalue(getvariable(m, :v))[1]
        cparcial[t] = getobjectivevalue(m)
    end


    writecsv("OUT_g1.csv",g1_vec)
    writecsv("OUT_g2.csv",g2_vec)
    writecsv("OUT_u.csv",u_vec)
    writecsv("OUT_v.csv",v_vec)
    writecsv("OUT_cparcial.csv",cparcial)
    # print(m_vec[1])
end

function backward!(prb::Problem, states::States, FCF, a_matrix::Array{Array{Float64,2},1})
    nStages = prb.nStages
    nBackwards = prb.nBackwards
    nForwards = prb.nForwards

    Cinf = 0
    for f in 1:nForwards
        for i in 0:(nStages - 2)
            t = nStages - i
            vec_γ1 = []
            vec_γ0 = []

            # calculates FCF approximations
            for b in 1:nBackwards
                a = a_matrix[t][b,:]
                dem = prb.dem[t][1]
                γ1, γ0 = subproblem(prb,states, FCF, a, dem, t, f, BACKWARD)

            push!(vec_γ1, γ1)
            push!(vec_γ0, γ0)
            end

            # updates FCF with expected value for approximations
            @show "back"
            avrg_γ1 = sum(1/nBackwards * vec_γ1)
            avrg_γ0 = sum(1/nBackwards * vec_γ0) - avrg_γ1 .* states.v[t,f,:]
            
            for pair in 1:length(avrg_γ1)
                push!(FCF[t, 1], (avrg_γ1[pair], avrg_γ0[pair]))
            end
            @show "backout"

            # 2 => first stage + 1 => FCF relative to first stage
            if t == 2
                # update Cinf
                Cinf = sum(1/nBackwards * vec_γ0)
            end
        end
    end

    return Cinf
end

function sddp()
    tolerance = 0.05

    # configuracao
    nStages = 25  # sempre um a mais
    nForwards = 10
    nBackwards = 10
    nHoursInStage = 1

    # dado de entrada
    # dem = 150 * ones(nStages)
    gcost = [100, 1500]
    umax = [300]
    vmax = [500]
    gmax = [200, 200]
    nGen = length(gmax)
    nHyd= length(umax)
    
    # Estimation of scenarios
    mat_dem = read_demand("hrload.csv")
    mat_inf = read_hinflw("inflow.csv")
    inf_forward, inf_backward = estimate(mat_inf, nForwards, nBackwards, nStages)
    writecsv("inflow_usado.csv",inf_forward)
    demand, nothing = estimate(mat_dem, 1, 1, nStages)
    writecsv("demanda_usada.csv",demand)
    
    # monta problema
    prb = Problem(nStages, nBackwards, nForwards, nHoursInStage, nGen, nHyd, demand, gmax, gcost, umax, vmax)
    inflw = Inflows(inf_forward, inf_backward)
    
    # states begin at stage 2, to avoid 0 index error
    states = States(prb.nStages, prb.nForwards, prb.nHyd)
    # initial storage
    states.v[1,:,:] = 50
    
    # FCF considers the FCF from last stage + 1
    nFCF = (prb.nStages + 1)
    FCF = [Tuple{Float64,Float64}[] for i in 1:nFCF, h in 1:nHyd]
    
    max_it = 20
    zsup_vec = []
    zinf_vec = []
    σ2 = 0
    Csup = 0
    Cinf = 0
    for i in 1:max_it
        Csup, σ2, Cinf = forward!(prb, states, FCF, inflw.forward)
        wrong = backward!(prb, states, FCF, inflw.backward)

        # build confidence 95% interval
        # IC_ub = Csup + Z * sqrt(σ2 / prb.nForwards)
        # IC_lb = Csup - Z * sqrt(σ2 / prb.nForwards)
        IC_ub = Csup + 2 * σ2
        IC_lb = Csup - 2 * σ2
        
        push!(zsup_vec, Csup)
        push!(zinf_vec, Cinf)

        # test for convergence
        if Cinf <= IC_ub &&  IC_lb <= Cinf && i > 3
            println("Solution found!")
            println("Zsup = $(Csup)")
            println("Zinf = $(Cinf)")
            break
        end
    end

    writecsv("OUT_zsup.csv",zsup_vec)
    writecsv("OUT_zinf.csv",zinf_vec)

    # simulate 
    simulation!(prb, states, FCF, inflw.forward)
end

sddp()