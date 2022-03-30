    subroutine Jumpterms(Calpha,ip,XI0,p1,U1,W1,cpts1,p2,U2,W2,cpts2,KN1,KN2,PN1,PN2,tol)
    implicit none
    !
    integer(kind=4)::ip,KN1,KN2,PN1,PN2,p1,p2
    real(kind=8)::Calpha,Xi0,tol
    real(kind=8),intent(in)::U1(0:KN1),W1(0:PN1),cpts1(0:PN1,1:2),U2(0:KN2),W2(0:PN2),cpts2(0:PN2,1:2)
    !
    integer(kind=4)::Jpoint(1:4),IEN(0:p1,1:2)
    real(kind=8)::pi
    real(kind=8)::XItemp(1:2),dRdXi(0:p1,1:2),dCi(1:2,1:2),JXi(1:2),tann(1:2,1:2),angle(1:2)
    logical::NOTJP
    logical::isJpoint(1:4)
    !
    pi = acos(-1.0D0)
    Jpoint(1) = 0
    Jpoint(2) = PN1
    Jpoint(3) = PN1+1
    Jpoint(4) = PN1+PN2+1
    isJpoint = (Jpoint == ip)
    NOTJP = .not. (any(isJpoint))
    !
    if ( NOTJP .and. (ip <= KN1) ) then
        Xitemp(1) = XI0 - tol*1.0D-2
        Xitemp(2) = Xi0 + tol*1.0D-2
        Call Ders1nrbasis(U1,KN1,Xitemp,2,W1,PN1,p1,dRdXi,IEN)
        Call RXicpts(dRdXi,p1,2,IEN,cpts1,PN1,2,dCi)
    elseif ( NOTJP .and. (ip > KN1) ) then
        Xitemp(1) = XI0 - tol*1.0D-2
        Xitemp(2) = Xi0 + tol*1.0D-2
        Call Ders1nrbasis(U2,KN2,Xitemp,2,W2,PN2,p2,dRdXi,IEN)
        Call RXicpts(dRdXi,p2,2,IEN,cpts2,PN2,2,dCi)
    elseif (isJpoint(1) .or. isJpoint(4)) then
        Xitemp(1) = U2(KN2) - tol*1.0D-2
        Xitemp(2) = U1(  0) + tol*1.0D-2
        Call Ders1nrbasis(U2,KN2,Xitemp(1:1),1,W2,PN2,p2,dRdXi(0:p2,1:1),IEN(0:p2,1:1))
        Call Ders1nrbasis(U1,KN1,Xitemp(2:2),1,W1,PN1,p1,dRdXi(0:p1,2:2),IEN(0:p1,2:2))
        Call RXicpts(dRdXi(0:p2,1:1),p2,1,IEN(0:p2,1:1),cpts2,PN2,2,dCi(1:1,1:2))
        Call RXicpts(dRdXi(0:p1,2:2),p1,1,IEN(0:p1,2:2),cpts1,PN1,2,dCi(2:2,1:2))
    elseif (isJpoint(2) .or. isJpoint(3)) then
        Xitemp(1) = U1(KN1) - tol*1.0D-2
        Xitemp(2) = U2(  0) + tol*1.0D-2
        Call Ders1nrbasis(U1,KN1,Xitemp(1:1),1,W1,PN1,p1,dRdXi(0:p1,1:1),IEN(0:p1,1:1))
        Call Ders1nrbasis(U2,KN2,Xitemp(2:2),1,W2,PN2,p2,dRdXi(0:p2,2:2),IEN(0:p2,2:2))
        Call RXicpts(dRdXi(0:p1,1:1),p1,1,IEN(0:p1,1:1),cpts1,PN1,2,dCi(1:1,1:2))
        Call RXicpts(dRdXi(0:p2,2:2),p2,1,IEN(0:p2,2:2),cpts2,PN2,2,dCi(2:2,1:2))
    else
        print*, "≈‰µ„ø’º‰Ω«≈–∂œ¥ÌŒÛ, Jumpterms.f90"
        print*, "Xi0 = ",XI0
        read(*,*)
        stop
    end if
    !
    JXi=sqrt(dCi(1:2,1)*dCi(1:2,1)+dCi(1:2,2)*dCi(1:2,2))
    tann(1:2,1) = dCi(1:2,1) / JXi
    tann(1:2,2) = dCi(1:2,2) / JXi
    !
    angle = atan2(tann(1:2,2),tann(1:2,1))

    Calpha = angle(2) - angle(1)
    if (Calpha >  pi) Calpha = Calpha - 2.0*pi
    if (Calpha < -pi) Calpha = Calpha + 2.0*pi
    Calpha = - (pi-Calpha) / (2.0*pi)

    
    if (Calpha < -1.0D0 .or. Calpha > 1.0D-10) then
        print*, "Calphaº∆À„¥ÌŒÛ"
        print*, "XI0 = ",Xi0
        read(*,*)
        stop
    end if 

    return
    end subroutine Jumpterms