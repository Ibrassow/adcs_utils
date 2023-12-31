using LinearAlgebra

include("transformations.jl")


function random_vector(μ::Vector, P::Matrix=I)
"""
    μ: mean vector
    P: covariance matrix
"""
    n = size(μ, 1)
    Δ = cholesky(P).L  

    return μ + randn(n,1) * Δ
end


function random_rotation_matrix(P::Matrix; Q0::Matrix=diagm(ones(3)))

    # Generate rnadom axis-angle vectors with desired covariance
    Δ = cholesky(P).L 
    ϕ = Δ * randn(3, 1)

    return Q0 * exp(skew_symmetric(ϕ))
end

function random_quaternion(P::Matrix)

    # Generate rnadom axis-angle vectors with desired covariance
    Δ = cholesky(P).L 
    ϕ = Δ * randn(3, 1)
    
    nϕ = norm(ϕ)
    return [ cos(nϕ/2) ; ϕ * sin(nϕ/2)/sqrt(ϕ'ϕ)]
end