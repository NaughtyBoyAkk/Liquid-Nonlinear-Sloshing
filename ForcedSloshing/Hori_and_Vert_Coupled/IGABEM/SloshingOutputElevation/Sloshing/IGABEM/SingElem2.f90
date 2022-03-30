    subroutine SingElem2(Ge,He,Ind,Xie,p,U,KN,W,WN,cpts,SXiN,ISing,tol)
    implicit none
    !  传递进来的变量声明
    integer(kind=4)::p,KN,WN,SXiN,ISing
    integer(kind=4)::Ind(0:p)
    real(kind=8)::tol
    real(kind=8)::Ge(0:p),He(0:p),Xie(1:2),U(0:KN),w(0:WN),cpts(0:WN,1:2)
    ! 其他变量声明
    integer(kind=4)::i,j,k,info_f,info_s
    integer(kind=4)::IENi(0:p,0:SXIN),JIn(1:SXIN+1),Ja(1:SXIN+1),Ia(1:SXIN+1),ipiv(1:SXiN)
    real(kind=8)::pi,DeltaSXi,Jfitemp,JCi0temp,ICi0temp
    real(kind=8)::Sxi(0:SXIN),RNi(0:p,0:SXIN),Ci(0:SXIN,1:2),dRNi(0:p,0:SXIN),dCi(0:SXIN,1:2),JXi(1:SXiN),outn(1:SxiN,1:2)
    real(kind=8)::deltax(1:SXiN),deltay(1:SXiN),deltar(1:SXiN),drdx(1:SxiN),drdy(1:SxiN),l12(1:SxiN,1:2),rl(1:SxiN)
    real(kind=8)::coffA(1:SXiN,1:SXiN),Ifitemp(1:SXiN),Hn(1:SXIN+1),Fn(1:SXIN+1),Jfi(0:p,1:SXiN),JCi0(0:p),Ifi(0:p,1:SXiN),ICi0(0:p)
    real(kind=8)::Jright(1:SXiN),Iright(1:SXiN),dJa(1:SXIN+1),dIa(1:SXIN+1)
    real(kind=8), external :: ddot
    !
    DeltaSXi=(Xie(2)-Xie(1))/dble(SXiN)
    if (ISing==1) then
        Sxi(0) = Xie(1) 
        do i=1,SXiN-1,1
            Sxi(i)=Sxi(i-1) + DeltaSXi
        end do
        Sxi(SXiN) = Xie(2)
    else
        Sxi(0) = Xie(2) 
        do i=1,SXiN-1,1
            Sxi(i)=Sxi(i-1) - DeltaSXi
        end do
        Sxi(SXiN) = Xie(1)
    end if
    !
    i=KN-p-1
    DeltaSXi = (Xie(1)+Xie(2)) / 2.0D0
    Call FindSpan(i, p, DeltaSXi, U, KN, j)
    do i=0,p,1
        Ind(i) = j-p+i +1
    end do
    do i=0,SXiN,1
        IENi(0:p,i) = Ind -1
    end do
    !
    Call Shell_nrbasis(U,KN,Sxi,SXiN+1,W,WN,p,RNi,IENi) 
    Call RXicpts(RNi,p,SXiN+1,IENi,cpts,WN,2,Ci)
    !
    if (ISing==1) Sxi(SxiN) = Sxi(SxiN) - tol * 1.0D-4
    Call Ders1nrbasis(U,KN,SXi,SXiN+1,W,WN,p,dRNi,IENi)
    Call RXicpts(dRNi,p,SXiN+1,IENi,cpts,WN,2,dCi)
    !!
    JXi=sqrt( dCi(1:SXiN,1)*dCi(1:SXiN,1) + dCi(1:SXiN,2)*dCi(1:SXiN,2) )
    outn(1:SXiN,1) =  dCi(1:SXiN,2) / JXi
    outn(1:SXiN,2) = -dCi(1:SXiN,1) / JXi
    deltax = Ci(1:SXiN,1) - Ci(0,1)
    deltay = Ci(1:SXiN,2) - Ci(0,2)
    deltar = sqrt(deltax*deltax + deltay*deltay)
    drdx = deltax/deltar
    drdy = deltay/deltar
    l12(1:SXiN,1) = -outn(1:SXiN,2)
    l12(1:SXiN,2) =  outn(1:SXiN,1) 
    rl = abs(drdx*l12(1:SXiN,1)+drdy*l12(1:SXiN,2))
    !!
    coffA(1:SXiN,1) = 1.0D0
    JIn(1) = 0
    JIn(SXiN+1) = SXiN
    do i = 2,SXiN,1
        coffA(1:SXiN,i) = coffA(1:SXiN,i-1) * deltar
        JIn(i) = i-1
    end do
    !!
    pi = acos(-1.0D0)
    Ja = 0 - JIn - 1
    Ia = JIn - 2 + 1
    dJa=dble(Ja)
    dIa=dble(Ia)
    !
    do k = 1,SXiN+1,1
        if (Ja(k)==0) then
            Hn(k) = log(deltar(SXiN))**2 / 2.0D0
        else
            Hn(k) = -(dJa(k) * log(deltar(SXiN))+1.0D0) / (dJa(k)*dJa(k)) / (deltar(SXiN)**dJa(k))
        end if
        if (Ia(k)==0) then
            Fn(k) = log(deltar(SXiN))
        else
            Fn(k) = deltar(SXiN)**dIa(k) / dIa(k)
        end if
    end do
    !
    !!
    Jfitemp  = 1.0D0 / (2.0D0*pi)
    Ifitemp  = ( deltax*outn(1:SXiN,1) + deltay*outn(1:SXiN,2) ) / (2.0D0 * pi)
    do i=1,SXiN,1
        Jfi(0:p,i) = Jfitemp    * RNi(0:p,i)
        Ifi(0:p,i) = Ifitemp(i) * RNi(0:p,i)
    end do
    JCi0temp = 1.0D0/(2.0D0*pi)
    JCi0 = JCi0temp * RNi(0:p,0)
    ICi0temp = 0.0D0
    ICi0 = ICi0temp * RNi(0:p,0)
    Call dgetrf(SXiN,SXiN,coffA,SXiN,ipiv,info_f)
    if (info_f /= 0) then
        print*, "奇异单元LU分解出错, SingElem2.f90"
        read(*,*)
        stop
    end if
    do i = 0, p,1
        Jright = (Jfi(i,1:SXiN)/rl-JCi0(i)) / deltar
        Call dgetrs('N',SXiN,1,coffA,SxiN,ipiv,Jright,SxiN,info_s)
        if (info_s /= 0) then
            print*, "奇异单元求解方程出错-J, SingElem2.f90"
            read(*,*)
            stop
        end if
        Ge(i) = ddot(SXIN,Jright,1,Hn(2:SXiN+1),1) + JCi0(i) * Hn(1)
        !
        Iright = (Ifi(i,1:SXiN)/rl-ICi0(i)) / deltar
        Call dgetrs('N',SXiN,1,coffA,SxiN,ipiv,Iright,SxiN,info_s)
        if (info_s /= 0) then
            print*, "奇异单元求解方程出错-I, SingElem2.f90"
            read(*,*)
            stop
        end if
        He(i) = ddot(SXIN,Iright,1,Fn(2:SXiN+1),1) + ICi0(i) * Fn(1)
    end do
    return
    end subroutine SingElem2