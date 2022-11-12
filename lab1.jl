using Plots
using LinearAlgebra

const fps = 24 # frames per second
const Δt = 1 / fps # s
const tmax = 30 # s
const times = 0:Δt:tmax # s
const L = 50# cm
const T0 = 293 # K
const kb = 1.38e-23 # J/K
const M = 1.63e-24 # kg
const N = 1000

function Maxwell2D(v, m=M, T=T0)
    return m * v / (kb * T) .* exp.(-m * v .^ 2 / (2 * kb * T))
end

function Maxwell3D(v, m=M, T=T0)
    return (m / (2 * π * kb * T))^(3 / 2) * 4π * v .^ 2 .* exp.(-m * v .^ 2 / (2 * kb * T))
end


function main()
    v0 = sqrt(3 * kb * T0 / M)
    r = L / 2 * ones(3, N)
    φ = rand(N) * 2π
    θ = rand(N) * π
    v = zeros(3, N)
    v[1, :] = v0 * cos.(φ) .* sin.(θ)
    v[2, :] = v0 * sin.(φ) .* sin.(θ)
    v[3, :] = v0 * cos.(θ)

    anim = Animation()
    v_hist = Animation()
    for time in times
        r = r + v * Δt
        too_left = r[1, :] .< 0
        too_right = r[1, :] .> L
        too_low = r[2, :] .< 0
        too_high = r[2, :] .> L
        too_near = r[3, :] .< 0
        too_far = r[3, :] .> L
        v[1, too_left] .= -v[1, too_left]
        v[2, too_low] .= -v[2, too_low]
        v[1, too_right] .= -v[1, too_right]
        v[2, too_high] .= -v[2, too_high]
        v[3, too_near] .= -v[3, too_near]
        v[3, too_far] .= -v[3, too_far]
        r[1, too_left] .= 0
        r[2, too_low] .= 0
        r[1, too_right] .= L
        r[2, too_high] .= L
        r[3, too_near] .= 0
        r[3, too_far] .= L

        for particle1 in 1:N
            for particle2 in particle1+1:N
                if particle1 == particle2
                    continue
                end
                if norm(r[:, particle1] - r[:, particle2]) < 0.7
                    v1 = v[:, particle1]
                    v2 = v[:, particle2]
                    r1 = r[:, particle1]
                    r2 = r[:, particle2]
                    v[:, particle1] = v1 - dot(v1 - v2, r1 - r2) / norm(r1 - r2)^2 * (r1 - r2)
                    v[:, particle2] = v2 - dot(v2 - v1, r2 - r1) / norm(r2 - r1)^2 * (r2 - r1)
                end
            end
        end

        vx = v[1, :]
        vy = v[2, :]
        vz = v[3, :]
        speed = sqrt.(vx .^ 2 .+ vy .^ 2 + vz .^ 2)
        histogram(speed, bins=100, legend=false, norm=true)
        # plot maxwell boltzmann distribution
        x = 0:0.1:2v0
        y = Maxwell3D(x)
        hist = plot!(x, y, lw=3, label="Maxwell-Boltzmann")
        frame(v_hist, hist)
        x = r[1, :]
        y = r[2, :]
        z = r[3, :]

        plt = scatter(x, y, z, xlims=(0, L), ylims=(0, L), zlims=(0, L), legend=false)
        frame(anim, plt)

    end
    gif(anim, "lab1.gif", fps=fps)
    gif(v_hist, "lab1_hist.gif", fps=fps)
end
main()