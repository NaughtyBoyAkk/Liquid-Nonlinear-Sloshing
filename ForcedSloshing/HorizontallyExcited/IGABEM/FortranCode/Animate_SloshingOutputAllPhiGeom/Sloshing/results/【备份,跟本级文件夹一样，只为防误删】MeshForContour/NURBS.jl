module NURBS
using SparseArrays
## #
function FindSpan(n::Int64,p::Int64,u::Float64,U::Vector{Float64})
    if u < U[1] || u > U[end]
        throw(DomainError(u, "u should be located in U"))
    end
    if u==U[n+2]
        knotSpanIndex= n
        return knotSpanIndex
    end
    low = p
    high = n+1
    mid = floor(Int,(low + high)/2)
    while u <U[mid+1] || u >= U[mid+2]
        if u < U[mid+1]
            high = mid
        else
            low = mid
        end
        mid = floor(Int,(low+high)/2)
    end
    knotSpanIndex = mid
    return knotSpanIndex
end
## #
function BasisFuns(i::Int64,u::Float64,p::Int64,U::Vector{Float64})
    Ni=zeros(Float64,p+1)
    left=copy(Ni)
    right=copy(Ni)
    Ni[1]=1.0
    for j=1:p
        left[j+1]=u-U[i+2-j]
        right[j+1]=U[i+j+1]-u
        saved=0.0
        for r=0:j-1
            temp=Ni[r+1]/(right[r+2]+left[j-r+1])
            Ni[r+1]=saved+right[r+2]*temp
            saved=left[j-r+1]*temp
        end
        Ni[j+1]=saved
    end
    return Ni
end
## #
function bsbasis(U::Vector{Float64},p::Int64,Xi::Vector{Float64})
    m=length(U)-1
    n=m-p-1
    k=length(Xi)
    N=spzeros(Float64,n+1,k)
    for i=1:k
        xi=Xi[i]
        j=FindSpan(n,p,xi,U)
        N[j-p+1:j+1,i]=BasisFuns(j,xi,p,U)
    end
    return N
end
## #
function nrbasis(U::Vector{Float64},w::Vector{Float64},Xi:: Vector{Float64})
    k=length(Xi)
    m=length(U)-1
    n=length(w)-1
    p=m-n-1
    N=bsbasis(U,p,Xi)
    #
    W=transpose(N)*w
    R=spzeros(Float64,n+1,k)
    for i=1:k
        R[:,i]=N[:,i].*w./W[i]
    end
    # R=N.*w./(transpose(w)*N) # 该公式破坏了稀疏性
    return R
end
## #
function DersBasisFuns(i::Int64,u::Float64,p::Int64,dn::Int64,U::Vector{Float64})
    if dn>p
        throw(DomainError(dn, "dn cannot be greater than p"))
    end
    ndu=zeros(Float64,p+1,p+1)
    left=zeros(Float64,p+1)
    right=copy(left)
    ders=zeros(Float64,dn+1,p+1)
    ndu[1,1]=1.0
    for j=1:p
        left[j+1]=u-U[i+2-j]
        right[j+1]=U[i+j+1]-u
        saved=0.0
        for r=0:j-1
            ndu[j+1,r+1]=right[r+2]+left[j-r+1]
            temp=ndu[r+1,j]/ndu[j+1,r+1]
            ndu[r+1,j+1]=saved+right[r+2]*temp
            saved=left[j-r+1]*temp
        end
        ndu[j+1,j+1]=saved
    end
    ders[1,:]=ndu[:,p+1]
    a=zeros(Float64,2,p+1)
    for r=0:p
        s1=0
        s2=1
        a[1,1]=1.0
        for k=1:dn
            d=0.0
            rk=r-k
            pk=p-k
            if r>=k
                a[s2+1,1]=a[s1+1,1]/ndu[pk+2,rk+1]
                d=a[s2+1,1]*ndu[rk+1,pk+1]
            end
            if rk>=-1
                j1=1
            else
                j1=-rk
            end
            if (r-1)<=pk
                j2=k-1
            else
                j2=p-r
            end
            for j=j1:j2
                a[s2+1,j+1]=(a[s1+1,j+1]-a[s1+1,j])/ndu[pk+2,rk+j+1]
                d+=a[s2+1,j+1]*ndu[rk+j+1,pk+1]
            end
            if r<=pk
                a[s2+1,k+1]=-a[s1+1,k]/ndu[pk+2,r+1]
                d+=a[s2+1,k+1]*ndu[r+1,pk+1]
            end
            ders[k+1,r+1]=d
            j=s1
            s1=s2
            s2=j
        end
    end
    r=p
    for k=1:dn
        for j=0:p
            ders[k+1,j+1]=ders[k+1,j+1]*r
        end
        r=r*(p-k)
    end
    return ders
end
## #
function Dersbsbasis(U::Vector{Float64},p::Int64,Xi::Vector{Float64},dn::Int64)
    m=length(U)-1
    n=m-p-1
    k=length(Xi)
    dNdξ=spzeros(Float64,n+1,k)
    for i=1:k
        xi=Xi[i]
        j=FindSpan(n,p,xi,U)
        dNdξ_loc=DersBasisFuns(j,xi,p,dn,U)
        dNdξ[j-p+1:j+1,i]=dNdξ_loc[dn+1,:]
    end
    return dNdξ
end
## #
function Dersnrbasis(U::Vector{Float64},w::Vector{Float64},
    Xi:: Vector{Float64},dn::Int64)
    m=length(U)-1
    n=length(w)-1
    p=m-n-1
    k=length(Xi)
    #
    N=bsbasis(U,p,Xi)
    dNdξ=Dersbsbasis(U,p,Xi,dn)
    #
    W=transpose(N)*w
    dWdξ=transpose(dNdξ)*w
    #
    dRdξ=spzeros(Float64,n+1,k)
    for i=1:k
        dRdξ[:,i]=(w.*(W[i].*dNdξ[:,i]-dWdξ[i].*N[:,i]))./(W[i]^2)
    end
    return dRdξ
end
## #
end
