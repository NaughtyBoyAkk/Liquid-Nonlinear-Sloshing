using Pkg
using Printf
using LinearAlgebra
using SparseArrays
using DelimitedFiles
using Plots
push!(LOAD_PATH, "./~")
include("NURBS.jl")
#
fU1 = open("knots1_p3.txt","r")
U1 = readdlm(fU1,Float64)
U1=[U1[1]; U1[:]; U1[end]]
close(fU1)
#
fU2 = open("knots2_p3.txt","r")
U2 = readdlm(fU2,Float64)
U2=[U2[1]; U2[:]; U2[end]]
close(fU2)
#
p = 3
#
m1 = length(U1) -1
m2 = length(U2) -1
#
PN1 = m1 - p - 1
PN2 = m2 - p - 1
w1 = ones(Float64,PN1+1)
w2 = ones(Float64,PN2+1)
#
plt1point=sort(unique([collect(LinRange(U1[1],U1[m1+1],250)); U1]))
plt2point=unique(U2)
R1=NURBS.nrbasis(U1,w1,plt1point)
R2=NURBS.nrbasis(U2,w2,plt2point)
#
tspan = collect(range(0.0,step=0.005,stop=12.0))
Nη = length(tspan)
η = zeros(Float64,Nη)
plt = plot(1,

           marker=4,
           xlim = (0.0, 12.0), ylim = (-0.1, 0.1),
           title = "Time history",
           legend = false  )
#
@gif for i=0:2399
    # i=0
    fc = open("cpts12_$(i).dat","r")
    cpts12 = readdlm(fc,Float64)
    close(fc)
    cpts1 = cpts12[1:PN1+1,1:2]
    cpts2 = cpts12[PN1+2:PN1+PN2+2,1:2]
    coors = [ transpose(R1)*cpts1;
              transpose(R2)*cpts2  ]
              #
    η[i+1] = coors[1,2]
    # 几何变化
    plot(coors[:,1], coors[:,2],marker=4,
    xlim = (-0.5, 0.5), ylim = (-0.65, 0.15), legend = false)
    #
end every 20
#
# 左侧角点波高
# @gif for i=0:2399
#     # i=0
#     fc = open("cpts12_$(i).dat","r")
#     cpts12 = readdlm(fc,Float64)
#     close(fc)
#     cpts1 = cpts12[1:PN1+1,1:2]
#     cpts2 = cpts12[PN1+2:PN1+PN2+2,1:2]
#     coors = [ transpose(R1)*cpts1;
#               transpose(R2)*cpts2  ]
#               #
#     η[i+1] = coors[1,2]
#     # 几何变化
#     # plot(coors[:,1], coors[:,2],marker=4,
#     # xlim = (-0.5, 0.5), ylim = (-0.65, 0.15), legend = false)
#     # 左侧角点波高
#     push!(plt,tspan[i+1],η[i+1])
# end every 10
