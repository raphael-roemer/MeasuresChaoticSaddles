using Plots
using LaTeXStrings
using Distributions

include("1-SystemDefinition.jl")

function mu_news(syst,rx,ry,ngridcells_xdir,sprinkle_n,N,m;
    tstep = 1)
    cx = range(rx[1],rx[2],ngridcells_xdir)
    cy = range(ry[1],ry[2],convert(Int64,ngridcells_xdir*((ry[2]-ry[1])/(rx[2]-rx[1]))))
    c_num = zeros(length(cx),length(cy))

    finals=zeros(sprinkle_n,2)
    rN = 0
    for i in 1:sprinkle_n
        #lock(lk)
        init = [rand(Uniform(rx[1],rx[2]),1)[1], rand(Uniform(ry[1],ry[2]),1)[1]]
        set_state!(syst, init)
        c = 0
        while ((rx[1]<=current_state(syst)[1]<=rx[2]) && (ry[1]<=current_state(syst)[2]<=ry[2]) && (c<=N))
            if c == m
                finals[i,1:2] = [current_state(syst)[1],current_state(syst)[2]]
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
            c_num[xx,yy] += 1
        end
        #unlock(lk)
    end

    c_mu = c_num/rN
    return c_mu,cx,cy
end


# set parameters
rx, ry = (0.2,1.9), (-0.408,0.408)
ngridcells_xdir*((ry[2]-ry[1])/(rx[2]-rx[1]))
rx_u, ry_u = (0.2,2.0), (-0.9,0.9)
tstep = 1


# measure on the saddle
ngridcells_xdir, sprinkle_n, N, m = 300, 300000000, 20, 2
c_mu, cx, cy = @time mu_news(syst,rx,ry,ngridcells_xdir, sprinkle_n,N,m)
Plots.heatmap(cx,cy,-log.(c_mu[:,:]'))




# # measure on saddle: inverse map

# ngridcells_xdir, ngridcells_zdir, sprinkle_n, N, m = 125, 125, 5000, 20, 20
# c_mui, cx, cy, cz = @time mu_news(syst_inv,rx,ry,rz,ngridcells_xdir, ngridcells_zdir, sprinkle_n,N,m)
# #Plots.heatmap(cx,cz,-log.(c_mu[5,:,:]'))
# c_mu_proj_xi = c_mui[1,:,:]
# c_mu_proj_yi = c_mui[:,1,:]
# for i in 1:ngridcells_zdir
#     c_mu_proj_xi[:,:] += c_mui[i,:,:]
#     c_mu_proj_yi[:,:] += c_mui[:,i,:]
# end
# Plots.heatmap(cx,cz,-log.(c_mu_proj_xi[:,:]'))
# Plots.heatmap(cy,cz,-log.(c_mu_proj_yi[:,:]'))



#fr = 4
#Plots.heatmap(cx[convert.(Int64,end/fr:2*end/3)],cy[convert.(Int64,end/fr:2*end/3)],log.(c_mu'[convert.(Int64,end/fr:2*end/3),convert.(Int64,end/fr:2*end/3)]))

 # measure on stable set
ngridcells_xdir_s, sprinkle_n_s, N_s, m_s = 200, 200000000, 16, 0
c_mu_s, cx_s, cy_s = @time mu_news(syst,rx,ry,ngridcells_xdir_s,sprinkle_n_s,N_s,m_s)
Plots.heatmap(cx_s,cy_s,log.(c_mu_s'))

# # measure on unstable set
ngridcells_xdir_u, sprinkle_n_u, N_u, m_u = 200, 200000000, 14, 14
c_mu_u, cx_u, cy_u = @time mu_news(syst,rx_u,ry_u,ngridcells_xdir_u,sprinkle_n_u,N_u,m_u)
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
#savefig("/Users/raphaelroemer/Desktop/measures.pdf")


# # plot all three measures in one plot
# c_mu_max = maximum(c_mu)
# c_mu_s_max = maximum(c_mu_s)
# c_mu_u_max = maximum(c_mu_u)

# p1 = Plots.heatmap(cx,cy,-log.(replace(x -> iszero(x) ? Inf : x, c_mu'./c_mu_max)), ylabel=L"y",xticks=false, xaxis=false ,c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match)
# p2 = Plots.heatmap(cx_s,cy_s,-log.(replace(x -> iszero(x) ? Inf : x, c_mu_s'./c_mu_s_max)), ylabel=L"y",xticks=false, xaxis=false,c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match)
# p3 = Plots.heatmap(cx_u,cy_u,-log.(replace(x -> iszero(x) ? Inf : x, c_mu_u'./c_mu_u_max)), ylabel=L"y", xlabel=L"x",c=reverse(cgrad(:inferno)),grid=false, leftmargin=:match)
# plot(p1, p2, p3, layout=(3,1), legend=false, size=(670,800))





