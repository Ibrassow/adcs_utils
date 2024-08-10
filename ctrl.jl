using LinearAlgebra
using PyCall 
pushfirst!(pyimport("sys")."path", ".")   

include("transformations.jl")
include("dynamics.jl")

function attitude_regulation(q_des, q, ω, Kp, Kd)
    # q_des: desired attitude
    # q: current attitude
    # ω: current angular velocity
    # Kp: proportional gain
    # Kd: derivative gain

    delta_q = L(conj(q)) * q_des 

    """if delta_q[1] == 0
        τ = - Kp * delta_q[2:end] - Kd * ω
    else
        τ = - sign(delta_q[1]) * Kp * delta_q[2:end] - Kd * ω
    end"""
    τ = -  Kp * delta_q[2:end] - Kd * ω

    return τ

end


# Sim test 

q_desired = [1.0, 0.0, 0.0, 0.0]
x0 = [0.434008,  0.685733,  -0.303806,  -0.499109, 0.0, 0.0, 0.0]

J = [0.1 0.0 0.0; 
    0.0 0.1 0.0; 
    0.0 0.0 0.1]
params = (
    J = J,
    invJ = inv(J)
    )

N = 1000

X = zeros(N, 7)
U = zeros(N-1, 3)
X[1, :] = x0
dt = 0.01
Kp = 1
Kd = 0.5

for i = 1:N-1
    U[i, :] = attitude_regulation(q_desired, X[i,1:4], X[i,5:7], Kp, Kd)
    X[i+1,:] = rk4mk(params, attitude_quat_rbd, X[i,:], U[i,:], dt)

end

plotter = pyimport("plotter")
importlib = pyimport("importlib")
importlib.reload(plotter)

plotter.plotAttitude(X, dt)


# MeshCat 

using MeshCat
using Colors 
using GeometryBasics
using CoordinateTransformations
using Rotations

vis = Visualizer()
setprop!(vis["/Background"], "top_color", Gray(0.0))
setprop!(vis["/Background"], "bottom_color", Gray(0.0))


c = 0.5
SC_box = HyperRectangle(Vec(-c/2, -c/2, -c/2), Vec(c, c, c))
setobject!(vis["SC"], SC_box)
dock_SC = HyperRectangle(Vec(-c/8, -c/8, -c/8) + [0.0, 0.0, c/2,], Vec(c/4, c/4, c/4))
setobject!(vis["SC"]["dock"], dock_SC, MeshPhongMaterial(color = RGBA(1,0,0,0.5)))

TT_SC = LinearMap(QuatRotation(X[1,1:4]))
settransform!(vis["SC"], TT_SC)



function update_attitude(vis1, q)
    settransform!(vis[vis1],  LinearMap(QuatRotation(q[1:4])))
end

sleep(5)

for i=1:N
    update_attitude("SC", X[i,:])
    sleep(dt)
end
