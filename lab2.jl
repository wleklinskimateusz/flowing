using Plots

function get_speed_profile(y::Vector{Float64})
    return 3 * sin.(π * (y) / (Δx * ny))
end
const Δx::Float64 = 1
const Δt::Float64 = 0.1
const tmax = 10
const D::Float64 = 4
const nx::Int64 = 100
const ny::Int64 = 100
const yrange::Vector{Float64} = 0:Δx:ny*Δx
const v::Vector{Float64} = get_speed_profile(yrange)


function get_sum(C::Matrix{Float64})
    return sum(sum(C))
end


function get_next_C_diffusion(C::Matrix{Float64}, i::Int64, j::Int64)::Float64
    if i == 1 && j == 1
        return C[i, j] + D * Δt * (C[i+1, j] + C[i, j+1] - 2 * C[i, j]) / Δx^2
    end
    if i == nx && j == ny
        return C[i, j] + D * Δt * (C[i-1, j] + C[i, j-1] - 2 * C[i, j]) / Δx^2
    end
    if i == 1 && j == ny
        return C[i, j] + D * Δt * (C[i+1, j] + C[i, j-1] - 2 * C[i, j]) / Δx^2
    end
    if i == nx && j == 1
        return C[i, j] + D * Δt * (C[i-1, j] + C[i, j+1] - 2 * C[i, j]) / Δx^2
    end
    if i == 1
        return C[i, j] + D * Δt / Δx^2 * (C[i+1, j] + C[i, j+1] + C[i, j-1] - 3 * C[i, j])
    end
    if j == 1
        return C[i, j] + D * Δt / Δx^2 * (C[i+1, j] + C[i-1, j] + C[i, j+1] - 3 * C[i, j])
    end

    if i == nx
        return C[i, j] + D * Δt / Δx^2 * (C[i-1, j] + C[i, j+1] + C[i, j-1] - 3 * C[i, j])
    end

    if j == ny
        return C[i, j] + D * Δt / Δx^2 * (C[i+1, j] + C[i-1, j] + C[i, j-1] - 3 * C[i, j])
    end

    return C[i, j] + D * Δt / (Δx^2) * (C[i+1, j] + C[i-1, j] + C[i, j+1] + C[i, j-1] - 4 * C[i, j])
end

function get_next_C_speed(C::Matrix{Float64}, i::Int64, j::Int64)::Float64
    if j == 1
        return C[i, j] + v[i] * Δt / Δx * (C[i, j+1] - C[i, j])
    end
    return C[i, j] + v[i] * Δt / Δx * (C[i, j-1] - C[i, j])
end




function main()
    C::Matrix{Float64} = zeros(nx, ny)
    C[10:90, 5:10] .= 1
    su = sum(sum(C))

    ani::Animation = Animation()
    times = 0:Δt:tmax
    S = zeros(length(times))
    for (ti, t) in enumerate(times)
        for i in 1:nx
            for j in 1:ny
                C[i, j] = get_next_C_diffusion(C, i, j)
                C[i, j] = get_next_C_speed(C, i, j)
            end
        end
        h = heatmap!(C, xlim=(0, 100), ylim=(0, 100), c=:jet, title="t = $t")
        frame(ani, h)
        sup = sum(sum(C))
        S[ti] = su - sup
        su = sup
    end
    println("saving")
    gif(ani, "lab2.gif", fps=1 / Δt)
    plot(times, S)
    savefig("S.png")
end
main()

# plot matrix
# heatmap(C, aspect_ratio=1, color=:viridis, legend=false)
# savefig("lab2.png")

