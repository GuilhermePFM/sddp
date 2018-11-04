struct Problem
    nStages::Int
    nBackwards::Int
    nForwards::Int
    nHours::Int
    nGen::Int
    nHyd::Int
    dem
    gmax
    gcost
    umax
    vmax
end

struct States
    v::Array{Float64, 3}
    States(nStages::Int, nFor, nHyd::Int) = new(zeros(nStages, nFor, nHyd))
end

struct Inflows
    forward::Array{Array{Float64,2},1}
    backward::Array{Array{Float64,2},1}
end

@enum Phase FORWARD=1 BACKWARD=2 SIMULATION=3