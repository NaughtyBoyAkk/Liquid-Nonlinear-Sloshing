    subroutine PartialTerms(knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2,nodes1,EN1,nodes2,EN2, &
        & Np,SXIN,etafp,phifp,keypoints,dphidx,dphidz,detadx)
    implicit none
    !
    integer(kind=4)::KN1,KN2,PN1,PN2,EN1,EN2,Np,SXIN,RK
    integer(kind=4)::keypoints(1:4)
    real(kind=8)::ti,analytical_eta
    real(kind=8)::knots1(0:KN1),knots2(0:KN2),cpts1(0:PN1,1:2),cpts2(0:PN2,1:2),weights1(0:PN1),weights2(0:PN2)
    real(kind=8)::nodes1(0:EN1),nodes2(0:EN2),etafp(0:EN1),phifp(0:EN1),dphidx(0:EN1),dphidz(0:EN1),detadx(0:EN1)
    !
    integer(kind=4)::p1,PNN
    integer(kind=4),allocatable::IENXi(:,:)
    real(kind=8),parameter::tol=1.0D-10
    real(kind=8),allocatable::phix(:),Veloy(:),etacpts(:,:),phicpts(:,:),RNg(:,:),dRgdXI(:,:),dCgdXi(:,:),JXI(:),outn(:,:)
    real(kind=8),allocatable::dphidXi(:),dphidn(:)
    !
    p1=KN1-PN1-1
    allocate(etacpts(0:PN1,1:1),phicpts(0:PN1,1:1))
    ! 计算液面波高的控制点
    Call CompuCpts(p1,knots1,Weights1,cpts1,KN1,PN1,Nodes1,EN1,Np,etafp,etacpts)
    ! 更新几何形状
    Call GeomUpdate(etacpts,cpts1,cpts2,PN1,PN2,keypoints)
    ! ! 边界条件-计算自由液面上phi边界条件的控制点
    Call CompuCpts(p1,knots1,Weights1,cpts1,KN1,PN1,Nodes1,EN1,Np,phifp,phicpts)
    !***********************************************************************************
    ! Laplace 2D 求解
    PNN=PN1+PN2+1 ! NN代表的总控制点数目是从0开始计的，所以这里是PN1+PN2+1，而不是PN1+PN2+2
    allocate(phix(0:PNN),Veloy(0:PNN))
    phix(0:PN1) = phicpts(0:PN1,1)
    phix(PN1+1:PNN) = 0.0D0
    Veloy = 0.0D0
    Call Laplace2D(knots1,knots2,Weights1,Weights2,cpts1,cpts2,nodes1,nodes2, & 
        & KN1,KN2,PN1,PN2,EN1,EN2,Np,SxiN,phix,Veloy,PNN)
    ! 计算偏导数
    allocate(outn(0:EN1,1:2),dphidXi(0:EN1),dphidn(0:EN1))
    allocate(RNg(0:p1,0:EN1),dRgdXI(0:p1,0:EN1),IENXi(0:p1,0:EN1),dCgdXi(0:EN1,1:2),JXI(0:EN1))
    !
    Call Ders1nrbasis(knots1,KN1,nodes1,EN1+1,Weights1,PN1,p1,dRgdXI,IENXI)
    Call RXicpts(dRgdXI,p1,EN1+1,IENXi,phix(0:PN1),PN1,1,dphidXi)
    !
    Call RXicpts(dRgdXI,p1,EN1+1,IENXi,cpts1,PN1,2,dCgdXi)
    JXi=sqrt(dCgdXi(0:EN1,1)*dCgdXi(0:EN1,1)+dCgdXi(0:EN1,2)*dCgdXi(0:EN1,2))
    outn(0:EN1,1) =  dCgdXi(0:EN1,2) / JXi
    outn(0:EN1,2) = -dCgdXi(0:EN1,1) / JXi
    Call nrbasis(knots1,KN1,nodes1,EN1+1,Weights1,PN1,p1,RNg,IENXI)
    Call RXicpts(RNg,p1,EN1+1,IENXi,Veloy(0:PN1),PN1,1,dphidn)
    !
    JXi = dCgdXi(0:EN1,1) * outn(0:EN1,2) - dCgdXi(0:EN1,2) * outn(0:EN1,1) !JXi 二次利用
    dphidx = ( outn(0:EN1,2) * dphidXi - dCgdXi(0:EN1,2) * dphidn ) / JXi
    dphidz = ( dCgdXi(0:EN1,1) * dphidn - outn(0:EN1,1) * dphidXI ) / JXi
    detadx = dCgdXi(0:EN1,2) / dCgdXi(0:EN1,1)
    !
    dphidx(0) = 0.0D0 
    dphidx(EN1) = 0.0D0 
    return
    end subroutine PartialTerms
