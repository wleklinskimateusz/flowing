using Plots
using LinearAlgebra

const Δt = 0.05 # s
const tmax = 7.5 # s
const times = 0:Δt:tmax # s
const L = 50# cm
const T = 293 # K
const kb = 1.38e-23 # J/K
const m = 1.63e-24 # kg
const N = 3000


function main()
    v0 = sqrt(3 * kb * T / m) # cm/s
    r = L / 2 * ones(2, N) #cm
    α = rand(N) * 2π # rad
    v = zeros(2, N) # cm/s
    v[1, :] = v0 * cos.(α)
    v[2, :] = v0 * sin.(α)
    anim = Animation()
    v_hist = Animation()
    for time in times
        r = r + v * Δt
        too_left = r[1, :] .< 0
        too_right = r[1, :] .> L
        too_low = r[2, :] .< 0
        too_high = r[2, :] .> L
        v[1, too_left] .= -v[1, too_left]
        v[2, too_low] .= -v[2, too_low]
        v[1, too_right] .= -v[1, too_right]
        v[2, too_high] .= -v[2, too_high]
        r[1, too_left] .= 0
        r[2, too_low] .= 0
        r[1, too_right] .= L
        r[2, too_high] .= L

        for particle1 in 1:N
            for particle2 in particle1+1:N
                if particle1 == particle2
                    continue
                end
                if norm(r[:, particle1] - r[:, particle2]) < 0.15
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
        speed = sqrt.(vx .^ 2 .+ vy .^ 2)
        histogram(speed, bins=100, legend=false, norm=true)
        # plot maxwell boltzmann distribution
        x = 0:0.1:2v0
        y = m * x / (kb * T) .* exp.(-m * x .^ 2 / (2 * kb * T))
        hist = plot!(x, y, lw=3, label="Maxwell-Boltzmann")
        frame(v_hist, hist)

        plt = scatter([r[1, :]], [r[2, :]], xlims=(0, L), ylims=(0, L), legend=false)
        frame(anim, plt)
    end
    gif(anim, "lab1.gif", fps=20)
    gif(v_hist, "lab1_hist.gif", fps=20)
end
main()