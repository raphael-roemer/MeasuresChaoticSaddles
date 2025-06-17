using DynamicalSystems


function expmap2(x, p, n)
    dx1 = p[1]*x[1]*exp(-(x[1])^2)
    dx2 = x[2] + p[2]*sin(x[2]) + p[3]*(x[1]^2-p[4])
    return SVector{2}(dx1,dx2)
end

function expmap2_inv(x, p, n)
    dx1 = p[1]*x[1]*exp(-(x[1])^2)
    dx2 = x[2] + p[2]*sin(x[2]) + p[3]*(x[1]^2-p[4])
    return SVector{2}(dx1,dx2)
end


function expmap2_deriv_for_lyap(x::Vector, p::Vector)
    dx1 = p[1]*exp(-(x[1])^2) - p[1]*x[1]*exp(-(x[1])^2)*(2*x[1])
    dx2 = 1 + p[2]*cos(x[2])
    return SVector{2}(dx1,dx2)
end

p=[4.2, 0.8, 0.4, 1.4]
ds4 = DiscreteDynamicalSystem(expmap2, zeros(2), p)
syst = ds4