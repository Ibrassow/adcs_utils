using LinearAlgebra
import Convex as cvx 
using Mosek
using MosekTools

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

landmarks = [rand(3) for i=1:100]



sz = size(landmarks, 1)

for i=1:sz
    landmarks[i] = landmarks[i] / norm(landmarks[i])
end
             
landmarks = Matrix(hcat(landmarks...)')


# Create respective unit vector measurements in body frame
measurements = []
for i=1:sz
    landmark = landmarks[i, :]
    v = true_att_Q' * (landmark - true_pos_eci)
    push!(measurements, v / norm(v))
end
measurements = Matrix(hcat(measurements...)')

Qcvx, _ = solve_wahba_cvx_relaxation(landmarks, measurements)
Qsvd = solve_wahba_svd(landmarks, measurements)
qdav, Qdav = solve_wahba_davenport_q_method(landmarks, measurements)