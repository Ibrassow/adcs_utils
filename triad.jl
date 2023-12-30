using LinearAlgebra

function TRIAD(n1::Vector, n2::Vector, b1::Vector, b2::Vector)
    # Classic TRIAD method
    # Need at least 2 unit vectors lnown in both inertial and body frames (measurement)
    # [n1 n2]' = Q [b1 b2]'

    # Normalization check
    n1 = n1 / norm(n1)
    n2 = n2 / norm(n2)
    b1 = b1 / norm(b1)
    b2 = b2 / norm(b2)


    # Inertial triad 
    t1 = n1
    t2 = cross(n1,n2) / norm(cross(n1,n2)) # Third lienarly independant vector 
    t3 = cross(t1,t2) / norm(cross(t1,t2))
    T = [t1 t2 t3]

    # Body triad
    w1 = b1
    w2 = cross(b1,b2) / norm(cross(b1,b2))
    w3 = cross(w1,w2) / norm(cross(w1,w2))
    W = [w1 w2 w3]

    # Determine attitude
    Q = T * W'
    
    return Q


end
