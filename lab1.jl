using Plots
using LinearAlgebra

const fps::Int64 = 110 # frames per second
const Δt::Float64 = 1 / fps # s
const tmax::Int64 = 2 # s
const times::Vector{Float64} = 0:Δt:tmax # s
const L::Float64 = 50# cm
const T0::Float64 = 293 # K
const kb::Float64 = 1.38e-23 # J/K
const M::Float64 = 1.63e-24 # kg
const N::Int64 = 2000
const vu::Float64 = 80

function Maxwell2D(v::Vector{Float64}, m::Float64=M, T::Float64=T0)::Vector{Float64}
    return m * v / (kb * T) .* exp.(-m * v .^ 2 / (2 * kb * T))
end

function Maxwell3D(v::Vector{Float64}, m=M, T=T0)::Vector{Float64}
    return (m / (2 * π * kb * T))^(3 / 2) * 4π * v .^ 2 .* exp.(-m * v .^ 2 / (2 * kb * T))
end

function calculate_particles_in_range(r::Matrix{Float64}, xmin::Float64, xmax::Float64)::Int64
    particles::BitVector = (r[1, :] .> xmin) .& (r[1, :] .< xmax)
    return sum(particles)
end


function main()
    v0::Float64 = sqrt(3 * kb * T0 / M)
    r::Matrix{Float64} = L / 2 * ones(3, N)
    φ::Vector{Float64} = rand(N) * 2π
    θ::Vector{Float64} = rand(N) * π
    v::Matrix{Float64} = zeros(3, N)
    v[1, :] = v0 * cos.(φ) .* sin.(θ)
    v[2, :] = v0 * sin.(φ) .* sin.(θ)
    v[3, :] = v0 * cos.(θ)

    anim = Animation()
    v_hist = Animation()
    x_min::Float64 = L
    x_max::Float64 = 1.1 * L
    particles_in_region = []
    for time in times
        bitmap = zeros(3, N)
        bitmap[1, :] .= 1
        r = r + (v + bitmap * vu) * Δt
        # too_left = r[1, :] .< 0
        # too_right = r[1, :] .> L
        too_low::BitVector = r[2, :] .< 0
        too_high::BitVector = r[2, :] .> L
        too_near::BitVector = r[3, :] .< 0
        too_far::BitVector = r[3, :] .> L
        # v[1, too_left] .= -v[1, too_left]
        v[2, too_low] .= -v[2, too_low]
        # v[1, too_right] .= -v[1, too_right]
        v[2, too_high] .= -v[2, too_high]
        v[3, too_near] .= -v[3, too_near]
        v[3, too_far] .= -v[3, too_far]
        # r[1, too_left] .= 0
        r[2, too_low] .= 0
        # r[1, too_right] .= L
        r[2, too_high] .= L
        r[3, too_near] .= 0
        r[3, too_far] .= L

        for particle1 in 1:N
            for particle2 in particle1+1:N
                if particle1 == particle2
                    continue
                end
                if norm(r[:, particle1] - r[:, particle2]) < 0.7
                    v1::Vector{Float64} = v[:, particle1]
                    v2::Vector{Float64} = v[:, particle2]
                    r1::Vector{Float64} = r[:, particle1]
                    r2::Vector{Float64} = r[:, particle2]
                    v[:, particle1] = v1 - dot(v1 - v2, r1 - r2) / norm(r1 - r2)^2 * (r1 - r2)
                    v[:, particle2] = v2 - dot(v2 - v1, r2 - r1) / norm(r2 - r1)^2 * (r2 - r1)
                end
            end
        end

        vx::Vector{Float64} = v[1, :]
        vy::Vector{Float64} = v[2, :]
        vz::Vector{Float64} = v[3, :]
        speed::Vector{Float64} = sqrt.(vx .^ 2 .+ vy .^ 2 + vz .^ 2)
        histogram(speed, bins=100, legend=false, norm=true)
        # plot maxwell boltzmann distribution
        x_data::Vector{Float64} = 0:0.1:2v0
        y = Maxwell3D(x_data)
        hist = plot!(x_data, y, lw=3, label="Maxwell-Boltzmann")
        frame(v_hist, hist)
        x::Vector{Float64} = r[1, :]
        y::Vector{Float64} = r[2, :]
        z::Vector{Float64} = r[3, :]

        plt = scatter(x, y, z, xlims=(-L, 5 * L), ylims=(0, L), zlims=(0, L), legend=false)
        frame(anim, plt)
        push!(particles_in_region, calculate_particles_in_range(r, x_min, x_max))

    end
    plot(times, particles_in_region, label="particles in region")
    savefig("lab1.png")
    gif(anim, "lab1.gif", fps=fps)
    gif(v_hist, "lab1_hist.gif", fps=fps)
end
main()