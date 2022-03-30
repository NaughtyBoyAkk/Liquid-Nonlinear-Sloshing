function GetXi(CB::Array{Float64,2},U::Vector{Float64},w::Vector{Float64},
    cpts::Array{Float64,2},tol::Float64=1.0e-8)
    #
    Nξ=size(CB,1)
    ξ0=zeros(Float64,Nξ)
    #
    Uu=unique(U)
    Uu[end] -= tol
    NU=length(Uu)
    RU=NURBS.nrbasis(U,w,Uu)
    CU=transpose(RU)*cpts
    #
    δx=repeat(CU[:,1],outer=(1,Nξ)).-transpose(CB[:,1])
    δy=repeat(CU[:,2],outer=(1,Nξ)).-transpose(CB[:,2])
    δr2=sqrt.(δx.^2+δy.^2)
    for i=1:Nξ
        minr,n1=findmin(δr2[:,i])
        if n1 <= length(Uu) && n1 >= 1 && abs(minr) < tol
            ξ0[i]=Uu[n1]
        elseif n1 <= length(Uu) && n1 >= 1 && abs(minr) > tol
            ξtemp=Uu[n1]+tol/2.0
            Rtemp=NURBS.nrbasis(U,w,[ξtemp])
            Ctemp=transpose(Rtemp)*cpts
            δtemp=sqrt((Ctemp[1]-CB[i,1])^2+(Ctemp[2]-CB[i,2])^2)
            if δr2[n1,i]>δtemp
                ξ12=Uu[n1:n1+1]
            else
                ξ12=Uu[n1-1:n1]
            end
            #
            while (ξ12[2]-ξ12[1])> tol
                R12=NURBS.nrbasis(U,w,ξ12)
                C12=transpose(R12)*cpts
                δx2 = (C12[:,1].-CB[i,1]).^2
                δy2 = (C12[:,2].-CB[i,2]).^2
                δr12 = sqrt.(δx2+δy2)
                if δr12[2]>δr12[1]
                    ξ12[2]=(ξ12[1]+ξ12[2])/2.0
                else
                    ξ12[1]=(ξ12[1]+ξ12[2])/2.0
                end
            end
            ξ0[i]=(ξ12[1]+ξ12[2])/2.0
        end
        Rξi = NURBS.nrbasis(U,w,ξ0[i:i])
        RC = transpose(Rξi) * cpts
        δv = norm((RC[:]-CB[i,:]),2)
        if δv > tol
            println((RC[:],CB[i,:]))
            throw("δv=$(δv),GetXi.jl未能找到正确解,请检查给定点是否属于给定几何")
        end
    end
    return ξ0
end
