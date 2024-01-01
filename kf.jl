using LinearAlgebra

mutable struct KalmanFilter

    
    x::Vector # Current estimate
    P::Matrix # Covariance matrix
    H::Matrix # Measurement matrix
    V::Matrix # Noise covariance matrix
    W::Matrix # Measurement covariance matrix
    A::Matrix # State transition matrix
    B::Matrix
    ctrl::Bool
    nx::Real
    nu::Real



    function KalmanFilter(x0::Vector, P0::Matrix, H::Matrix, V::Matrix, W::Matrix, A::Matrix, B::Matrix=nothing)
        nx = size(x0, 1)
        if B !== nothing
            nu = size(B, 1)
        else
            ctrl=false
        end
        # TODO Check all sizes
        new(x0, P0, HC, V, W, A, B, ctrl, nx, nu)
    end

end

function predict(kf::KalmanFilter, u::Vector=zeros(kf.nu))
    if ctrl == false
        kf.x = kf.A * kf.x
    else
        kf.x = kf.A * kf.x + kf.B * u
    end


    kf.P = kf.A * kf.P * kf.A' + kf.W

    return kf.x, kf.P
end

function update(kf::KalmanFilter, y_mes::Vector)

    # Innovation + Covariance
    z = y_mes - kf.H * kf.x
    S = kf.H * kf.P * kf.H' + kf.V

    # Kalman gain
    K = kf.P * kg.H' * inv(S)

    # update
    kf.x = kf.x + K * z
    kf.P = (I - K * kf.H) * kf.P * (I - K * kf.H)' + K * V * K'

    return kf.x, kf.P

end