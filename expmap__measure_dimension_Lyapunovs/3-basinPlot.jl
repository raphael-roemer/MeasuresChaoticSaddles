using Plots
using LaTeXStrings
using Distributions
using ProgressMeter

#include("1-SystemDefinition.jl")

function basin_plot_vals(syst,rx,ry,sprinkle_n,N;
    tstep = 1)
    inits=zeros(5,sprinkle_n) # x, y, attractor_number, x_final, y_final
    inits[1,:] = rand(Uniform(rx[1],rx[2]),sprinkle_n)
    inits[2,:] = rand(Uniform(ry[1],ry[2]),sprinkle_n)

    tstep_total=tstep*N
    #Threads.@threads 
    for i in 1:sprinkle_n
        set_state!(syst, inits[1:2,i])
        # step forwards exactly tstep_total
        step!(syst, tstep_total) # add ,true) when using a continous system
        inits[4:5,i] = current_state(syst)
        if inits[5,i] > 1.3
            inits[3,i] = 1
        else
            inits[3,i] = 0
        end
    end
    return inits
end


# set parameters and prepare plot
sprinkle_n_bas, N_bas = 22000, 50
rx_bas = (0.01,2.0)
ry_bas = (-4,4.0)
tstep = 1
plot_vals = @time basin_plot_vals(syst,rx_bas,ry_bas,sprinkle_n_bas,N_bas, tstep=tstep)

#choose colorscheme
basincolors_pal = palette([palette(:sun)[end-1], :lightblue], 2);
attracolors_pal = palette([:darkblue, :darkblue], 2);
#attracolors_pal = palette([palette(:sun)[end-2], :darkblue], 2)
basincolors=basincolors_pal[convert.(Int8,plot_vals[3,:]).+1];
fra = 1
attracolors=attracolors_pal[convert.(Int8,plot_vals[3,1:convert(Int64,sprinkle_n_bas/fra)]).+1]

#Plots.plot(plot_vals[4,:], plot_vals[5,:],ms=1,alpha=0.99,markercolor=basincolors,markerstrokecolor=basincolors,seriestype=:scatter,label=false,grid=false) 
Plots.plot(plot_vals[1,:], plot_vals[2,:],ms=1.5,alpha=0.39,markercolor=basincolors,markerstrokecolor=basincolors,seriestype=:scatter,label=false,grid=false,size=(900,600)) 
Plots.plot!(plot_vals[4,1:convert(Int64,sprinkle_n_bas/fra)],plot_vals[5,1:convert(Int64,sprinkle_n_bas/fra)],ms=1.5,alpha=0.1,markercolor=attracolors,markerstrokecolor=attracolors,seriestype=:scatter,ylim=(-4.5,4.5), label=false)
Plots.plot!(plot_vals[4,1:2],plot_vals[5,1:2],ms=1.5,alpha=1,markercolor=:blue,markerstrokecolor=:blue,seriestype=:scatter,ylim=(-4.5,4.5), label=L"attractors")


# plot saddle and stable set in the plot
saddle = zeros(2,length(c_mu[1,:])*length(c_mu[:,1]))
stable = zeros(2,length(c_mu[1,:])*length(c_mu[:,1]))

countstable = 1
countsaddle = 1


for i in 1:length(c_mu[:,1])
    for j in 1:length(c_mu[1,:])
        if !iszero(c_mu[i,j])
            saddle[:,countsaddle] = [cx[i],cy[j]]
            countsaddle += 1
        end
        if !iszero(c_mu_s[i,j])
            stable[:,countstable] = [cx[i],cy[j]]
            countstable += 1
        end

    end
end
Plots.plot!(stable[1,1:countstable],stable[2,1:countstable], seriestype=:scatter,ms=1.5,alpha=0.39,markercolor=:lime,markerstrokecolor=:lime,label=false,grid=false, xlabel=L"x", ylabel=L"y")
Plots.plot!(stable[1,1:2],stable[2,1:2], seriestype=:scatter,ms=1.5,alpha=1,markercolor=:lime,markerstrokecolor=:lime,label=L"stable \ set",grid=false, xlabel=L"x", ylabel=L"y")
Plots.plot!(saddle[1,1:countsaddle],saddle[2,1:countsaddle], seriestype=:scatter,ms=1.5,alpha=0.39,markercolor=:green,markerstrokecolor=:green,label=false,grid=false)
Plots.plot!(saddle[1,1:2],saddle[2,1:2], seriestype=:scatter,ms=1.5,alpha=1,markercolor=:green,markerstrokecolor=:green,label=L"chaotic \ saddle",grid=false)

savefig("/Users/raphaelroemer/Desktop/basins.pdf")