struct Problem
    nStages::Int
    nBackwards::Int
    nForwards::Int
end

struct States
    v::Array{Float64, 1}
    States(nStages::Int) = new(zeros(nStages))
end

@enum Phase FORWARD=1 BACKWARD=2 SIMULATION=3