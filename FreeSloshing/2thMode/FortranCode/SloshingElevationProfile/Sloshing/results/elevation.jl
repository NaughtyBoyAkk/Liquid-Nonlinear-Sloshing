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
plt1point=sort(unique([collect(LinRange(U1[1],U1[m1+1],1000)); U1]))
R1=NURBS.nrbasis(U1,w1,plt1point)
#
for i in [532; 551; 570; 589; 608]
    fc = open("cpts1_$(i).dat","r")
    cpts12 = readdlm(fc,Float64)
    close(fc)
    cpts1 = cpts12[1:PN1+1,1:2]
    coors = transpose(R1)*cpts1;
    fid = open("elevation$(i).dat","w")
    writedlm(fid,coors)
    close(fid)
    #
end
