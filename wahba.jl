using LinearAlgebra
import Convex as cvx 
using Mosek
using MosekTools
import ForwardDiff as FD
using Printf

include("transformations.jl")



function solve_wahba_cvx_relaxation(landmarks::Matrix, measurements::Matrix)

    ################ Variable ################

    Q = cvx.Variable(3,3)
    sz = size(landmarks, 1)


    ################ Objective function ################

    B = zeros(3,3)

    # Building attitude profile matrix
    # assume equal weight for now
    for k=1:sz
        B += 1.0 * measurements[k,:] * landmarks[k,:]'
    end

    cost = -tr(B*Q)
    prob = cvx.minimize(cost)

    ################ Constraints ################
    id = diagm(ones(3))
    prob.constraints +=  [[id Q'; Q id] == cvx.Semidefinite(6)]

    ################ Output ################

    cvx.solve!(prob, Mosek.Optimizer; silent_solver = false)

    if (prob.status != cvx.MathOptInterface.OPTIMAL) && (prob.status != cvx.MathOptInterface.ALMOST_OPTIMAL)
        error("Convex.jl problem failed to solve.")
    end

    Q = Matrix(Q.value)
    final_cost = prob.optval

    return Q, final_cost

end

function solve_wahba_svd(landmarks::Matrix, measurements::Matrix)

    sz = size(landmarks, 1)
    B = zeros(3,3)

    # Building attitude profile matrix
    # assume equal weight for now
    for k=1:sz
        B += 1.0 * measurements[k,:] * landmarks[k,:]'
    end

    sv = svd(B)

    md = diagm(ones(3))
    md[3,3] = det(sv.U) * det(sv.V)
    nQb = sv.U * md * sv.Vt
    return nQb

end

function solve_wahba_davenport_q_method(landmarks::Matrix, measurements::Matrix)

    sz = size(landmarks, 1)
    B = zeros(3,3)
    z = zeros(3)

    # Building attitude profile matrix
    # assume equal weight for now
    for k=1:sz
        B += 1.0 * measurements[k,:] * landmarks[k,:]'
        z += 1.0 * (cross(landmarks[k,:], measurements[k,:]))
    end


    # Construct Davenport K matrix
    #K = [tr(B) z'; z (B + B' - tr(B)*I)]
    K = [(B + B' - tr(B)*I) z; z' tr(B) ]

    eigval, eigvec = eigen(K)
    _, ind = findmax(eigval)

    qq = eigvec[ind, :]/norm(eigvec[ind, :])
    q = zeros(4)
    q[1] = qq[4]
    q[2:4] = qq[1:3]
    q = q / norm(q)

    #q = eigvec[ind, :]/norm(eigvec[ind, :])
    nQb = dcm_from_q(q)

    return q, nQb

end

function angular_error(q1,q2)
    q = qconj(q1) ⊙ q2
    return (2 * atan(norm(q[2:4]), q[1]))
end


function solve_wahba_multiplicative_gauss_newton(landmarks::Matrix, measurements::Matrix, q0::Vector; tol::Real=1e-6, max_iter::Real=1000)

    sz = size(landmarks, 1)
    

    function residual(q)
        r = zeros(eltype(q),sz*3) # loss fct is r'r
        for i=1:sz
            r[ (i-1)*3 + 1 : (i-1)*3 + 3] = measurements[i,:] - dcm_from_q(q) * landmarks[i, :]
        end
        return r
    end

    function loss(q)
        r = residual(q)
        return r' * r
    end

    # Jacobian w.r.t quaternion
    Jquat(qval) = FD.ForwardDiff.jacobian(q -> residual(q), qval)
    # Attitude Jacobian
    G(q) = L(q) * [zeros(eltype(q), 1, 3); Matrix{eltype(q)}(I, 3, 3)]

    # Initial guess
    q = q0

    for iter = 1:max_iter

        jac = Jquat(q) * G(q)

        # Three-parameter update step
        ϕ = - inv(jac' * jac) * jac' * residual(q)

        curr_loss = loss(q)

        α = 1.0
        new_loss = 0.0

        # Line search 
        for l=1:20

            # line search step
            qn = L(q) * quat_from_mrp(α * ϕ)

            # new loss
            new_loss = loss(qn)

            if new_loss < curr_loss
                # sufficient decrease ~ armijo? 
                q = copy(qn)
                break
            else
                # decrease step length
                α = 0.5 * α
            end


        end

        ΔL = new_loss - curr_loss
        if (iter % 10) == 0
            @printf "iter       new_loss            ΔL            α            \n"
            @printf "-------------------------------------------------\n"
        end

        @printf("%3d    %10.3e      %10.3e      %9.2e       \n", iter, new_loss, ΔL, α)

        if new_loss < tol
            break
        end

    end

    return q
end

# Ground truth
true_pos_eci = [0, 0.0, 0.0]
true_att_q = [0.9725208, 0.0864657, 0.1296985, -0.1729313]
true_att_Q = dcm_from_q(true_att_q)
println("Ground truth position\\n", true_pos_eci)
println("Ground truth attitude\\n", true_att_Q)
println("###################################################")

# Landmarks position in ECI
landmarks = [[2.5, -5.75, 10],
             [-1.5, -4.5, 8],
             [0.3, -5.0, 5],
             [2.0, -3.75, 4],
             [1.0, -4.2, 3],
             [1.0, 0.0, 6]]

landmarks = [rand(3) for i=1:10]



sz = size(landmarks, 1)

for i=1:sz
    landmarks[i] = landmarks[i] / norm(landmarks[i])
end
             
landmarks = Matrix(hcat(landmarks...)')


# Create respective unit vector measurements in body frame
measurements = []
for i=1:sz
    landmark = landmarks[i, :]
    v = true_att_Q' * (landmark - true_pos_eci) #+ randn(3)*0.1
    push!(measurements, v / norm(v))
end
measurements = Matrix(hcat(measurements...)')

Qcvx, _ = solve_wahba_cvx_relaxation(landmarks, measurements)
errQcvx = Matrix(true_att_Q' * Qcvx)
err_cvx_deg = (180 / π) *  norm(unskew_symmetric(log(errQcvx)))
@show err_cvx_deg

Qsvd = solve_wahba_svd(landmarks, measurements)
errQsvd = Matrix(true_att_Q' * Qsvd)
err_svd_deg = (180 / π) *  norm(unskew_symmetric(log(errQsvd)))
@show err_svd_deg

qdav, Qdav = solve_wahba_davenport_q_method(landmarks, measurements)
errQdav = Matrix(true_att_Q' * Qdav)
err_dav_deg = (180 / π) *  norm(unskew_symmetric(log(errQdav)))
@show err_dav_deg

qgn = solve_wahba_multiplicative_gauss_newton(landmarks, measurements, qdav)
Qgn = dcm_from_q(qgn)
errQgn = Matrix(true_att_Q' * Qgn)
err_gn_deg = (180 / π) *  norm(unskew_symmetric(log(errQgn)))
@show err_gn_deg