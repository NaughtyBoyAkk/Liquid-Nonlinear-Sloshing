    subroutine NonSingElem2(Ge,He,Ind,Np,C0,Xie,p,U,KN,w,WN,cpts)
    implicit none
    ! 传递进来的变量声明
    integer(kind=4)::Np,p,KN,WN
    integer(kind=4)::Ind(0:p)
    real(kind=8)::Ge(0:p),He(0:p),C0(1:2),Xie(1:2),U(0:KN),w(0:WN),cpts(WN,2)
    ! 其他变量声明
    integer(kind=4)::i
    integer(kind=4)::IEN(0:p,1:Np)
    real(kind=8)::Jeta,pi
    real(kind=8)::xp(1:Np),wf(1:Np),Xigauss(1:Np),RNg(0:p,1:Np),Cgauss(1:Np,2),dRNg(0:p,1:Np),dCgauss(1:Np,2)
    real(kind=8)::JXi(1:Np),outn(1:Np,2),deltax(1:Np),deltay(1:Np),deltar2(1:Np),g(1:Np),h(1:Np)
    pi=acos(-1.0D0)
    Call GaussPointsWeightFactors(-1.0D0,1.0D0,xp,wf,Np)
    Jeta = (Xie(2)-Xie(1))/2.0D0
    Xigauss = (1.0D0 - xp) / 2.0D0 * Xie(1) + (1.0D0 + xp) / 2.0D0 * Xie(2)
    Call nrbasis(U,KN,Xigauss,Np,W,WN,p,RNg,IEN)
    if (IEN(p,Np)-IEN(p,1) /= 0) then
        print*, "请检查奇异单元内Gauss点的分布, NonSingElem2.f90"
        read(*,*)
        stop
    end if
    Call RXicpts(RNg,p,Np,IEN,cpts,WN,2,Cgauss)
    !
    Call Ders1nrbasis(U,KN,Xigauss,Np,W,WN,p,dRNg,IEN)
    if (IEN(p,Np)-IEN(p,1) /= 0) then
        print*, "请检查奇异单元内Gauss点的分布, NonSingElem2.f90"
        read(*,*)
        stop
    end if
    Call RXicpts(dRNg,p,Np,IEN,cpts,WN,2,dCgauss)
    JXi=sqrt(dCgauss(1:Np,1)*dCgauss(1:Np,1)+dCgauss(1:Np,2)*dCgauss(1:Np,2))
    outn(1:Np,1) =  dCgauss(1:Np,2)/JXi
    outn(1:Np,2) = -dCgauss(1:Np,1)/JXi
    deltax = Cgauss(1:Np,1) - C0(1)
    deltay = Cgauss(1:Np,2) - C0(2)
    deltar2 = deltax*deltax + deltay*deltay
    g = log(sqrt(deltar2))/(2.0D0*pi)
    h = (deltax*outn(1:Np,1) + deltay*outn(1:Np,2))/deltar2/(2.0D0*pi)
    Ge = 0.0D0
    He = 0.0D0
    do i=1,Np,1
        Ge = Ge + RNg(0:p,i)*g(i)*JXi(i)*Jeta*wf(i)
        He = He + RNg(0:p,i)*h(i)*JXi(i)*Jeta*wf(i)
    end do
    Ind = IEN(0:p,Np) + 1
    return
    end subroutine NonSingElem2