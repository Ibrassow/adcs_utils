using LinearAlgebra
using PyCall

if !isdefined(Main, :plotter) 
    plotter = pyimport("plotter")
    importlib = pyimport("importlib")
end

importlib.reload(plotter)

include("dynamics.jl")



params = (
    J = [5000.0 0 0; 0 9000 0; 0 0 12000],
)

q0 = [1.0, 0, 0, 0]
ω0 = [1.1300, -0.2300, 0.8]


N = 100
x = zeros(eltype(q0), N, 7)
x[1, :] = [q0 ; ω0]


for i=1:N-1
    x[i+1, :] = rk4mk(params, attitude_quat_rbd, x[i, :], u0, dt)
end

plotter.plot3axis(x[:,5:7], dt, true)


## We want to spin abour the y-axis 
Jeff = 1.2 * params.J[3,3]
params.J[2,2] = Jeff

## Use choose ρ such that Jeff > J_33, e.g. 1.2 * J_33