    subroutine CompuCpts(p,U,W,cpts,KN,PN,Nodes,EN,Np,etafp,etacpts)
    implicit none
    !
    integer(kind=4)::p,KN,PN,EN,Np
    real(kind=8)::U(0:KN),W(0:PN),cpts(0:PN,1:2),Nodes(0:EN),etafp(0:EN),etacpts(0:PN,1:1)
    !
    integer(kind=4)::i,j,k,Ngauss
    integer(kind=4)::IENi(0:p),ipiv(1:PN+1)
    real(kind=8)::xp0(1:Np),wf0(1:Np),RTR(1:PN+1,1:PN+1),F(1:PN+1)
    integer(kind=4),allocatable::IENXi(:,:)
    real(kind=8),allocatable::XIgauss(:),etagauss(:),wf(:),Jeta(:),RNg(:,:),dRgdXI(:,:),dCgdXi(:,:),Jxi(:)
    !
    Call GaussPointsWeightFactors(-1.0D0,1.0D0,xp0,wf0,Np)
    !
    Ngauss = Np*EN
    allocate(XIgauss(1:Ngauss),etagauss(1:Ngauss),wf(1:Ngauss),Jeta(1:Ngauss),IENXi(0:p,1:Ngauss))
    do i=1,EN,1
        j = (i-1) * Np + 1
        k = i * Np
        wf(j:k) = wf0(1:Np)
        Jeta(j:k) = ( nodes(i) - nodes(i-1) ) / 2.0D0
        XIgauss(j:k) = (1.0D0 - xp0) / 2.0D0 * nodes(i-1) + (1.0D0 + xp0) / 2.0D0 * nodes(i)
        etagauss(j:k)= (1.0D0 - xp0) / 2.0D0 * etafp(i-1) + (1.0D0 + xp0) / 2.0D0 * etafp(i)
    end do
    !
    allocate(RNg(0:p,1:Ngauss),dRgdXI(0:p,1:Ngauss),dCgdXi(1:Ngauss,1:2),Jxi(1:ngauss))
    Call Ders1nrbasis(U,KN,Xigauss,Ngauss,W,PN,p,dRgdXI,IENXI)
    Call RXicpts(dRgdXI,p,Ngauss,IENXi,cpts,PN,2,dCgdXi)
    JXi=sqrt(dCgdXi(1:ngauss,1)*dCgdXi(1:ngauss,1)+dCgdXi(1:ngauss,2)*dCgdXi(1:ngauss,2))
    Call nrbasis(U,KN,Xigauss,Ngauss,W,PN,p,RNg,IENXI)
    !
    RTR = 0.0D0
    F = 0.0D0
    do i=1,Ngauss,1
        IENi = IENXi(0:p,i) + 1
        do j=0,p,1
            k = IENXI(j, i) + 1
            RTR(IENi, k) = RTR(IENi, k) + RNg(0:p,i) * RNg(j,i) * Jeta(i) * JXi(i) * wf(i)
        end do
        F(IENi) = F(IENi) + RNg(0:p,i) * etagauss(i) * Jeta(i) * JXi(i) * wf(i)
    end do
    Call dgetrf(PN+1,PN+1,RTR,PN+1,ipiv,i)
    if (i /= 0) then
        print*, "边界条件求解，LU分解出错, BCs1.f90"
        read(*,*)
        stop
    end if
    Call dgetrs('N',PN+1,1,RTR,PN+1,ipiv,F,PN+1,j)
    if (j /= 0) then
        print*, "边界条件求解，求解出错, BCs1.f90"
        read(*,*)
        stop
    end if
    etacpts(0:PN,1) = F(1:PN+1)
    return
    end subroutine CompuCpts