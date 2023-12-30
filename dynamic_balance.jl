using LinearAlgebra

include("dynamics.jl")


params = (
    J = [100.0 0 0; 0 120 0; 0 0 85],
)

q0 = [1.0, 0, 0, 0]
ω0 = [0.01, 0.02, -0.3]
x0 = [q0 ; ω0]

x = zeros(eltype(x0), N, 7)

Jeff = 1.2 * params.J