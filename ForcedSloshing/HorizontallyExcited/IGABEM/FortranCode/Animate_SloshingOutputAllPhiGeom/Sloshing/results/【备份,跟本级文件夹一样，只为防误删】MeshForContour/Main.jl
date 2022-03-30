using Pkg
using Printf
using LinearAlgebra
using SparseArrays
using DelimitedFiles
using Plots
push!(LOAD_PATH, "./~")
include("GetXi.jl")
include("NURBS.jl")
include("genGeom.jl")
include("Post.jl")
include("gauss.jl")
## 画云图用的FE网格
L = 0.9
H = 0.6
#
Ex = 120
Ey = 80
EN = Ex *Ey
#
Nx = Ex + 1
Ny = Ey + 1
NN = Nx * Ny
#
Nodes = zeros(Float64,NN,2)
Elems = zeros(Int32,EN,4)
#
Nodes[1:Ny,1] .= -L / 2.0
Nodes[1:Ny,2] = collect(range(-H, 0.0, length = Ny))
for i = 1:Ex
    nn = collect(i*Ny+1 : (i+1)*Ny)
    Nodes[nn,1] .= L / Ex * i - L / 2.0
    Nodes[nn,2] = Nodes[1:Ny,2]
    for j = 1:Ey
        k = (i-1) * Ey + j
        n1 = (i-1) * Ny + j
        n2 = n1 + Ny
        n3 = n2 + 1
        n4 = n1 +1
        Elems[k,1:4] = [n1; n2; n3; n4]
    end
end
# Plots.plot(Nodes[:,1],Nodes[:,2],
# legend=:none,
# color=:black,
# marker=4)
## 找到云图所用的FE网格中，自由液面对应的节点编号和xy坐标ηxy
ηNodesNum = collect(Nx*Ny:-Ny:Ny)
ηxy = Nodes[ηNodesNum,1:2]
# 导入IGABEM的初始构型
p = 3
k1 = "knots1_p$(p).txt"
c1 = "cpts1_p$(p).txt"
k2 = "knots2_p$(p).txt"
c2 = "cpts2_p$(p).txt"
U1,cpts1,w1,U2,cpts2,w2 = genGeom(k1,c1,k2,c2)
## 求解FE网格中自由液面对应的IGABEM初始构型的ξ值
ηξ = GetXi(ηxy,U1,w1,cpts1)
# 各时刻构型
m1 = length(U1) -1
m2 = length(U2) -1
#
PN1 = m1 - p - 1
PN2 = m2 - p - 1
## 导入变形后的IGABEM构型，并求解波高η
deformedNodes = copy(Nodes)
i = "8.5000"
fc = open("cpts12_$(i).dat","r")
cpts12 = readdlm(fc,Float64)
close(fc)
ηfcpts = cpts12[1:PN1+1,1:2]
wallcpts = cpts12[PN1+2:end,1:2]
Rη=NURBS.nrbasis(U1,w1,ηξ)
η = transpose(Rη) * ηfcpts
## 根据波高η和FE网格，将FE网格进行变形，更新FE各节点的纵坐标
for j=1:Nx
    nn = collect((j-1)*Ny+1 : j*Ny)
    α = (Nodes[nn,2] .- Nodes[nn[1],2]) ./ (Nodes[nn[end],2]-Nodes[nn[1],2])
    Δy = α .* ( η[Nx-j+1,2] - Nodes[nn[end],2])
    deformedNodes[nn,2] = Nodes[nn,2] .+ Δy
end
## 画出节点查看
# Plots.plot(deformedNodes[:,1],deformedNodes[:,2],
# legend=:none,
# color=:black,
# marker=4)
## 将网格节点信息输出文件备用
# fid = open("Nodes$(i).dat","w")
# writedlm(fid,deformedNodes)
# close(fid)
# fid = open("Elems$(i).dat","w")
# writedlm(fid,Elems)
# close(fid)
## 开始求解FE中各节点处的速度势
# 找到FE网格中容器壁对应的节点
tankwall = [collect(Ny:-1:2); collect(1:Ny:Ny*Ex+1); collect(Ny*Ex+2:Ny*Nx)]
# 边界节点，只是在调试时使用了一次
# bNNum = [ηNodesNum; tankwall]
# 自由液面的节点编号为ηNodesNum，容器壁的节点编号为tankwall
## 接下来读入IGABEM速度势和法向速度的控制点
fid = open("phicpts$(i).dat","r")
phicpts = readdlm(fid, Float64)
close(fid)
phicpts1 = phicpts[1:PN1+1]
phicpts2 = phicpts[PN1+2:end]
#
fid = open("vcpts$(i).dat","r")
vcpts = readdlm(fid, Float64)
close(fid)
vcpts1 = vcpts[1:PN1+1]
vcpts2 = vcpts[PN1+2:end]
## 开始FE节点循环
# 特别注意，变形后的IGABEM构型的控制点为ηfcpts 和 wallcpts
Np=5
xp,wf=gauss(Float64,Np)
NN1 = length(ηNodesNum)
NN2 = length(tankwall)
phi = zeros(Float64,NN)
for j=1:NN
    # j = 176
    coors = deformedNodes[j:j,1:2]
    in1 = findfirst((ηNodesNum .- j).== 0)
    in2 = findfirst((tankwall .- j).== 0)
    if in1==nothing && in2==nothing #域内点
        Gp1,Hp1 = Post(xp,wf,coors,U1,w1,ηfcpts)
        Gp2,Hp2 = Post(xp,wf,coors,U2,w2,wallcpts)
        phi[j:j] = (Hp1*phicpts1+Hp2*phicpts2) - (Gp1*vcpts1+Gp2*vcpts2)
    elseif in1==nothing && in2 !==nothing #位于容器壁
        ξB = GetXi(coors,U2,w2,wallcpts)
        RξB= NURBS.nrbasis(U2,w2,ξB)
        phi[j:j] = transpose(RξB) * phicpts2
    elseif in1 !==nothing && in2==nothing #位于自由表面
        ξB = GetXi(coors,U1,w1,ηfcpts)
        RξB= NURBS.nrbasis(U1,w1,ξB)
        phi[j:j] = transpose(RξB) * phicpts1
    else # 自由液面与容器壁的交点处，用来判断速度势是否连续
        ξB = GetXi(coors,U2,w2,wallcpts)
        RξB= NURBS.nrbasis(U2,w2,ξB)
        phi2 = transpose(RξB) * phicpts2
        #
        ξB = GetXi(coors,U1,w1,ηfcpts)
        RξB= NURBS.nrbasis(U1,w1,ξB)
        phi1 = transpose(RξB) * phicpts1
        if (norm(phi1-phi2,2) > 1.0e-6)
            throw("节点速度势不连续，请检查原因")
        else
            phi[j:j] = transpose(RξB) * phicpts1
        end
    end
end
res=open("phiIGAContour$(i).dat","w")
@printf(res,"TITLE = \"Contour\"\n")
@printf(res,"VARIABLES = \"x\",\"z\",\"<greek>f</greek>\"\n")
@printf(res,"ZONE DATAPACKING = POINT, NODES = %8d, ELEMENTS = %8d, ZONETYPE = FEQUADRILATERAL\n",NN,EN)
writedlm(res,[deformedNodes[:,1] deformedNodes[:,2] phi])
writedlm(res,Elems)
close(res)
