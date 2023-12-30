using LinearAlgebra
using PyCall

if !isdefined(Main, :plotter) 
    plotter = pyimport("plotter")
    importlib = pyimport("importlib")
end

importlib.reload(plotter)


params = (
    J = [10000.0 0 0; 0 9000 0; 0 0 12000],
)

q0 = [1.0, 0, 0, 0]
ω0 = [1.1300, -0.2300, 0.8]


N = 1000
x = zeros(eltype(q0), N, 7)
x[1, :] = [q0 ; ω0]

u0 = zeros(3)
dt = 0.1

for i=1:N-1
    x[i+1, :] = rk4mk(params, attitude_quat_rbd, x[i, :], u0, dt)
end

plotter.plot3D(x[:,5:7], center=true, equal_axis=true)

