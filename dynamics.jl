using LinearAlgebra
using ExponentialAction

include("transformations.jl")


function quat_kinematics(q::Vector, ω::Vector; H::Matrix=[zeros(eltype(q), 1, 3); Matrix{eltype(q)}(I, 3, 3)])
    """
    q: (w,x,y,z)
    """
    q̇ = 0.5 * L(q) * H * ω
    return q̇
end

function attitude_quat_rbd(params::NamedTuple, x::Vector, τ::Vector)
    """
    Rotational dynamics only
    """
    q = x[1:4] # Attitude (quaternions)
    ω = x[5:7] # Angular Velocity
    q̇ = quat_kinematics(q, ω)
   
    ωdot = params.J \ (τ -  skew_symmetric(ω) * params.J * ω)

    return [q̇ ; ωdot]
end


function euler_rotational_dynamics(params::NamedTuple, ω::Vector, τ::Vector)
    return params.J \ (τ - cross(ω, params.J * ω))
end

function euler_rotational_dynamics_inv(params::NamedTuple, ω::Vector, τ::Vector)
    return params.invJ * (τ - cross(ω, params.J * ω))
end


function rk4(params::NamedTuple, dynamics::Function, x::Vector, u::Vector, dt::Real)
    k1 = dt * dynamics(params, x, u)
    k2 = dt * dynamics(params, x + k1 * 0.5, u)
    k3 = dt * dynamics(params, x + k2 * 0.5, u)
    k4 = dt * dynamics(params, x + k3, u)
    x = x + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return x
end



function rk4mk(P::NamedTuple, dynamics::Function, xk::Vector, uk::Vector, dt::Real)
    H = [zeros(eltype(xk), 1, 3); Matrix{eltype(xk)}(I, 3, 3)]

    k1 = dt * dynamics(P, xk, uk)
    k2 = dt * dynamics(P, xk + k1 * 0.5, uk)
    k3 = dt * dynamics(P, xk + k2 * 0.5, uk)
    k4 = dt * dynamics(P, xk + k3, uk)

    f = (xk + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4))

    q = xk[1:4]
    ω = xk[5:7]

    ## Using ExponentialAction - ForwardDiff friendly https://github.com/sethaxen/ExponentialAction.jl 
    f[1:4] = expv(1, R( H * (0.5 * dt * ω + ((dt)/6)*(k1[5:7] + k2[5:7] + k3[5:7]))), q) 

    return f
end


