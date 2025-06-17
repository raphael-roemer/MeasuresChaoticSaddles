using DynamicalSystems
using NonlinearSolve
using StaticArrays
using Plots

function expmap3_bakers(x, p, n)
    if x[1] <= p[1]
        dx1 = x[1]/p[1]
        dx2 = x[2]*p[5]
        dx3 = x[3] + p[2]*sin(x[3]) + p[3]*(x[1]^2-p[4])
    else
        dx1 = (x[1]-p[1])/(1 - p[1])
        dx2 = (1 - p[5])*x[2]+p[5]
        dx3 = x[3] + p[2]*sin(x[3]) + p[3]*(x[1]^2-p[4])
    end
    return SVector{3}(dx1,dx2,dx3)
end

f(u, m) = u.+ m[1].*sin.(u).+m[2].*(m[4]^2 .-m[3]).-m[5]
u0 = @SVector[0.0]

function expmap3_bakers_inv(x, p, n)
    if x[2] < p[5]
        dx1 = x[1]*p[1]
        dx2 = x[2]/p[5]
        m = [p[2], p[3], p[4],dx1, x[3]]
        probN = NonlinearProblem(f, u0, m)
        dx3 = solve(probN, NewtonRaphson(), reltol = 1e-12, abstol=1e-15)[1]
    else
        dx1 = x[1]*(1 - p[1]) + p[1]
        dx2 = (x[2]-p[5])/(1 - p[5])
        m = [p[2], p[3], p[4],dx1, x[3]]
        probN = NonlinearProblem(f, u0, m)
        dx3 = solve(probN, NewtonRaphson(), reltol = 1e-12, abstol=1e-15)[1]
    end
    return SVector{3}(dx1,dx2,dx3)
end

#x = [0.1,0.1,0.8]
#expmap3_bakers_inv(expmap3_bakers(x,p,0),p,0)

p=[0.7, 0.4, 0.4, 0.5, 0.51] # asymmetric
#p=[0.5, 0.4, 0.4, 0.5, 0.5] # symmetric
ds1 = DiscreteDynamicalSystem(expmap3_bakers, zeros(3), p)
syst = ds1

dsi1 = DiscreteDynamicalSystem(expmap3_bakers_inv, zeros(3), p)
syst_inv = dsi1



