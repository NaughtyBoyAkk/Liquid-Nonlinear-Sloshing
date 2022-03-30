function Post(ηg::Vector{Float64},wf::Vector{Float64},C0::Array{Float64,2},
    U::Vector{Float64},w::Vector{Float64},cpts::Array{Float64,2})
    #
    # 生成单元的knots表
    nodes=unique(U)
    NN=length(nodes)
    elems=[nodes[1:NN-1] nodes[2:NN]]
    wf=repeat(wf,outer=NN-1)
    # 计算ξ到η等参变换的雅可比矩阵的行列式
    Np=length(ηg)
    Jη=repeat((elems[:,2]-elems[:,1])./2.0,inner=Np)
    totalElems=Vector(1:1:NN-1)
    gaussnumber=reshape(Vector(1:1:(NN-1)*Np),Np,NN-1)
    # 等参变换，将正则空间中的Gauss points转换到参数空间
    # 存为一个列向量，每Np个为一个单元的所有gauss points
    # 所以ξg长度为Np*(NN-1)
    ξg=reshape([(1.0 .-ηg)./2.0 (1.0 .+ηg)./2.0]*transpose(elems),:)
    # 计算各Gauss points对应的物理坐标Cg
    # 也就是各单元的gauss points对应的物理位置
    Rg=NURBS.nrbasis(U,w,ξg)
    Cg=transpose(Rg)*cpts
    # 计算曲线在物理空间的点r到ξ等参变换的雅可比矩阵的行列式
    dRgdξ=NURBS.Dersnrbasis(U,w,ξg,1)
    dCgdξ=transpose(dRgdξ)*cpts
    Jξ=sqrt.(dCgdξ[:,1].^2+dCgdξ[:,2].^2)
    #
    nξg=[dCgdξ[:,2] -dCgdξ[:,1]]./Jξ
    #
    m=length(U)-1
    n=length(w)-1
    p=m-n-1
    #
    Nξ0=size(C0,1)
    Gp=zeros(Float64,Nξ0,n+1)
    Hp=zeros(Float64,Nξ0,n+1)
    #
    for i=1:Nξ0
        δx=Cg[:,1].-C0[i,1]
        δy=Cg[:,2].-C0[i,2]
        δr2=δx.^2+δy.^2
        # 需要注意的是，由于本程序利用稀疏矩阵存储的全局基底矩阵R
        # 所以不需要像传统BEM或FEM那样组装，直接相加即可
        green=log.(sqrt.(δr2))./(2.0*pi)
        Gp[i,:]=sum(transpose(Rg).*green.*wf.*Jξ.*Jη,dims=1)
        #
        dgreendn=sum([δx δy].*nξg,dims=2)./(δr2.*(2.0*pi))
        Hp[i,:]=sum(transpose(Rg).*dgreendn.*wf.*Jξ.*Jη,dims=1)
        #
    end
    return Gp,Hp
end
