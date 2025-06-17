using Plots
using LaTeXStrings
using Distributions
#using ProgressMeter
#using Base.Threads

#include("1-SystemDefinition.jl")




# nthreads()
# Threads.threadid()
#Threads.@threads for i in 1:6
#    local te = zeros(nthreads())
#    print(Threads.threadid())
#end
#rN = zeros(1,nthreads())

# idea: ask how many threads, and then create an array with this dimension for 
# init, c, c_num, rN
# aggregate the values from the different kernels

function mu_news(syst,rx,ry,rz,ngridcells_xdir,ngridcells_zdir,sprinkle_n,N,m;
    tstep = 1)
    cx = range(rx[1],rx[2],ngridcells_xdir)
    cy = range(ry[1],ry[2],convert(Int64,ngridcells_xdir*((ry[2]-ry[1])/(rx[2]-rx[1]))))
    cz = range(rz[1],rz[2],ngridcells_zdir)
    c_num = zeros(length(cx),length(cy),length(cz))

    #lk = ReentrantLock()
    finals=zeros(sprinkle_n,3)
    rN = 0
    for i in 1:sprinkle_n
        #lock(lk)
        init = [rand(Uniform(rx[1],rx[2]),1)[1], rand(Uniform(ry[1],ry[2]),1)[1], rand(Uniform(rz[1],rz[2]),1)[1]]
        set_state!(syst, init)
        c = 0
        while ((rx[1]<=current_state(syst)[1]<=rx[2]) && (ry[1]<=current_state(syst)[2]<=ry[2]) && (rz[1]<=current_state(syst)[3]<=rz[2]) && (c<=N))
            if c == m
                finals[i,1:3] = [current_state(syst)[1],current_state(syst)[2], current_state(syst)[3]]
            end
            # step forwards exactly tstep
            step!(syst, tstep) # add ,true) when using a continous system
            c += 1 
        end
        if c-1 == N
            rN += 1
            # if not escaped then find where it was at step m
            xx=convert(Int64,floor((length(cx)-1)*(finals[i,1]-rx[1])/(rx[2]-rx[1]))+1)
            yy=convert(Int64,floor((length(cy)-1)*(finals[i,2]-ry[1])/(ry[2]-ry[1]))+1)
            zz=convert(Int64,floor((length(cz)-1)*(finals[i,3]-rz[1])/(rz[2]-rz[1]))+1)
            c_num[xx,yy,zz] += 1
        end
    end

    c_mu = c_num/rN
    return c_mu,cx,cy,cz
end


# set parameters
rx, ry, rz = (0.0,1.0), (0.0,1.0), (-1.0,1.0)

# measure on saddle

ngridcells_xdir, ngridcells_zdir, sprinkle_n, N, m = 125, 125, 70000000, 30, 0;
c_mu, cx, cy, cz = @time mu_news(syst,rx,ry,rz,ngridcells_xdir, ngridcells_zdir, sprinkle_n,N,m)
#Plots.heatmap(cx,cz,-log.(c_mu[5,:,:]'))
c_mu_proj_x = c_mu[1,:,:];
c_mu_proj_y = c_mu[:,1,:];
for i in 1:ngridcells_zdir
    c_mu_proj_x[:,:] += c_mu[i,:,:]
    c_mu_proj_y[:,:] += c_mu[:,i,:]
end;
#Plots.heatmap(cx,cz,-log.(c_mu_proj_x[:,:]'))
Plots.heatmap(cy,cz,-log.(c_mu_proj_y[:,:]'))




# measure of the inverse map

ngridcells_xdir, ngridcells_zdir, sprinkle_n, N, m = 125, 125, 5000, 30, 30;
c_mui, cx, cy, cz = @time mu_news(syst_inv,rx,ry,rz,ngridcells_xdir, ngridcells_zdir, sprinkle_n,N,m)
#Plots.heatmap(cx,cz,-log.(c_mu[5,:,:]'))
c_mu_proj_xi = c_mui[1,:,:];
c_mu_proj_yi = c_mui[:,1,:];
for i in 1:ngridcells_zdir
    c_mu_proj_xi[:,:] += c_mui[i,:,:]
    c_mu_proj_yi[:,:] += c_mui[:,i,:]
end;
#Plots.heatmap(cx,cz,-log.(c_mu_proj_xi[:,:]'))
Plots.heatmap(cy,cz,-log.(c_mu_proj_yi[:,:]'))


#fr = 4
#Plots.heatmap(cx[convert.(Int64,end/fr:2*end/3)],cy[convert.(Int64,end/fr:2*end/3)],log.(c_mu'[convert.(Int64,end/fr:2*end/3),convert.(Int64,end/fr:2*end/3)]))

# measure on stable set in forward time
ngridcells_xdir, ngridcells_zdir, sprinkle_n, N, m = 125, 125, 90000000, 30, 0;
c_mu_sForward, cx, cy, cz = @time mu_news(syst,rx,ry,rz,ngridcells_xdir, ngridcells_zdir, sprinkle_n,N,m)
c_mu_sForward_proj_x = c_mu_sForward[1,:,:];
c_mu_sForward_proj_y = c_mu_sForward[:,1,:];
for i in 1:ngridcells_zdir
    c_mu_sForward_proj_x[:,:] += c_mu_sForward[i,:,:]
    c_mu_sForward_proj_y[:,:] += c_mu_sForward[:,i,:]
end;
#Plots.heatmap(cx,cz,-log.(c_mu_proj_x[:,:]'))
Plots.heatmap(cy,cz,-log.(c_mu_sForward_proj_y[:,:]'),xlabel="y",ylabel="z")
#savefig("/Users/raphaelroemer/Desktop/SymmStabForward.pdf")
savefig("/Users/raphaelroemer/Desktop/ASymmStabForward.pdf")

# measure on unstable set in inverse time

ngridcells_xdir, ngridcells_zdir, sprinkle_n, N, m = 125, 125, 5000, 30, 30;
c_mui_uBackward, cx, cy, cz = @time mu_news(syst_inv,rx,ry,rz,ngridcells_xdir, ngridcells_zdir, sprinkle_n,N,m)
#Plots.heatmap(cx,cz,-log.(c_mu[5,:,:]'))
c_mui_uBackward_proj_x = c_mui_uBackward[1,:,:];
c_mui_uBackward_proj_y = c_mui_uBackward[:,1,:];
for i in 1:ngridcells_zdir
    c_mui_uBackward_proj_x[:,:] += c_mui_uBackward[i,:,:]
    c_mui_uBackward_proj_y[:,:] += c_mui_uBackward[:,i,:]
end;
#Plots.heatmap(cx,cz,-log.(c_mu_proj_xi[:,:]'))
Plots.heatmap(cy,cz,-log.(c_mui_uBackward_proj_y[:,:]'),xlabel="y",ylabel="z")
#savefig("/Users/raphaelroemer/Desktop/SymmUnsBackward.pdf")
savefig("/Users/raphaelroemer/Desktop/ASymmUnsBackward.pdf")






ngridcells_xdir_s, sprinkle_n_s, N_s, m_s = 300, 20000000, 16, 0
c_mu_s, cx_s, cy_s = @time mu_new(syst,rx,ry,ngridcells_xdir_s,sprinkle_n_s,N_s,m_s)
Plots.heatmap(cx_s,cy_s,log.(c_mu_s'))

# # measure on unstable set
ngridcells_xdir_u, sprinkle_n_u, N_u, m_u = 300, 20000000, 14, 14
c_mu_u, cx_u, cy_u = @time mu_new(syst,rx_u,ry_u,ngridcells_xdir_u,sprinkle_n_u,N_u,m_u)
Plots.heatmap(cx_u,cy_u,log.(c_mu_u'))


# plot all three measures in one plot
c_mu_max = maximum(-log.(replace(x -> iszero(x) ? Inf : x, c_mu')))
c_mu_s_max = maximum(-log.(replace(x -> iszero(x) ? Inf : x, c_mu_s')))
c_mu_u_max = maximum(-log.(replace(x -> iszero(x) ? Inf : x, c_mu_u')))

p1 = Plots.heatmap(cx,cy,-log.(replace(x -> iszero(x) ? Inf : x, c_mu'))./c_mu_max, ylabel=L"y",xticks=false, xaxis=false ,c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match, size=(2700,900))
#savefig("/Users/raphaelroemer/Desktop/measures1.pdf")
p2 = Plots.heatmap(cx_s,cy_s,-log.(replace(x -> iszero(x) ? Inf : x, c_mu_s'))./c_mu_s_max, ylabel=L"y",xticks=false, xaxis=false,c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match,size=(2700,900))
#savefig("/Users/raphaelroemer/Desktop/measures2.pdf")
p3 = Plots.heatmap(cx_u,cy_u,-log.(replace(x -> iszero(x) ? Inf : x, c_mu_u'))./c_mu_u_max, ylabel=L"y", xlabel=L"x",c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match,size=(2700,900))
#savefig("/Users/raphaelroemer/Desktop/measures3.pdf")
Plots.plot(p1, p2, p3, layout=(3,1), legend=false,size=(900,1000)) #size=(900,1000),size=(670,800))
#savefig("/Users/raphaelroemer/Desktop/measures.png")


# # plot all three measures in one plot
# c_mu_max = maximum(c_mu)
# c_mu_s_max = maximum(c_mu_s)
# c_mu_u_max = maximum(c_mu_u)

# p1 = Plots.heatmap(cx,cy,-log.(replace(x -> iszero(x) ? Inf : x, c_mu'./c_mu_max)), ylabel=L"y",xticks=false, xaxis=false ,c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match)
# p2 = Plots.heatmap(cx_s,cy_s,-log.(replace(x -> iszero(x) ? Inf : x, c_mu_s'./c_mu_s_max)), ylabel=L"y",xticks=false, xaxis=false,c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match)
# p3 = Plots.heatmap(cx_u,cy_u,-log.(replace(x -> iszero(x) ? Inf : x, c_mu_u'./c_mu_u_max)), ylabel=L"y", xlabel=L"x",c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match)
# plot(p1, p2, p3, layout=(3,1), legend=false, size=(670,800))





