using DataFrames
using Distributions #Random

function read_demand(file)
    df = readtable(file,skipstart=1)

    return Matrix(df[[:BPA_IDMT]])
end

function read_hinflw(file)
    df = readtable(file)

    return Matrix(df[[:Inflow]])
end

# inflows structure: array of matrices
# array dimension: stages
# matrix dimensions: line is scenario and column is agent
function estimate(mat, nScenFor, nScenBack, nStages)
    nAgent = size(mat)[2]
    inflow_forw = [zeros(nScenFor, nAgent) for t in 1:nStages]
    inflow_back = [zeros(nScenFor+nScenBack, nAgent) for t in 1:nStages]

    for s in 1:nScenFor

        for t in 1:nStages
            idx = nStages*(s-1) + t
            inflow_forw[t][s,:] = mat[idx,:]
            inflow_back[t][s,:] = mat[idx,:]
        end
    end

    # fit noise
    vec = mat[1:744,:]
    σ =  std(vec)
    μ = mean(vec)
    d = LogNormal(μ, σ)
    
    for sback in 1:nScenBack
        for t in 1:nStages
            idx = nStages*(sback-1) + t
            ϵ = rand(d)
            inflow_back[t][nScenFor+sback,:] = log(ϵ) #+ mat[idx,:]
        end
    end
    return inflow_forw, inflow_back
end

# function estimation(mat, nSeries)
#     nStages, nHyd = size(mat)
#     sint_series = ones(nStages, nSeries, nHyd)

#     sint_series = vec .* sint_series

#     # fit noise
#     σ =  std(vec)
#     μ = mean(vec)
#     d = Normal(0.0, σ/maximum(vec))
#     ϵ = rand(d, length(vec), nSeries)

#     return mat2array(sint_series + ϵ)
# end

# function mat2array(mat)
#     lines, cols =size(mat)

#     arr = []
#     for i in 1:lines
#         arr_temp = []
#         for j in 1:cols
#             push!(arr_temp, mat[i,j])
#         end
#         push!(arr, arr_temp)
#     end
#     return arr
# end


# mat_dem = read_demand("hrload.csv")
# mat_inf = read_hinflw("inflow.csv")

# estimate(mat_inf, 3, 3, 10)
# estimate(mat_dem, 3, 3,744)