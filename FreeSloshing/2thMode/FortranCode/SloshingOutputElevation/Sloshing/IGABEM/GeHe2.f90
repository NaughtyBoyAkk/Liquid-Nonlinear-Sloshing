    subroutine GeHe2(Ge,He,Ind,Xi0,C0,Xie,p2,knots2,KN2,weights2,PN2,cpts2,l,Np,SXiN,tol)
    implicit none
    !���ݽ����ı�������
    integer(kind=4)::KN2,PN2,p2,Np,SXiN
    integer(kind=4)::Ind(0:p2)
    real(kind=8)::Xi0,tol
    real(kind=8)::knots2(0:KN2),cpts2(0:PN2,1:2),weights2(0:PN2)
    real(kind=8)::Ge(0:p2),He(0:p2),Xie(1:2),C0(1:2)
    logical::l(1:10)
    ! ������������
    integer(kind=4)::ISing
    integer(kind=4)::Ind2(0:p2)
    real(kind=8)::Xie2(1:2),Ge2(0:p2),He2(0:p2)
    !
    if (count(l) /= 1) then
        print*, "����������жϴ���, GeHe2.f90"
        read(*,*)
        stop
    elseif (l( 2)) then ! �ڶ��������ڵĵ�һ�ڵ�����
        ISing=1
        Call SingElem2(Ge,He,Ind,Xie,p2,knots2,KN2,weights2,PN2,cpts2,SXiN,ISing,tol)
    elseif (l( 4)) then ! �ڶ��������ڵĵڶ��ڵ�����
        ISing=2
        Call SingElem2(Ge,He,Ind,Xie,p2,knots2,KN2,weights2,PN2,cpts2,SXiN,ISing,tol)
    elseif (l( 6)) then ! �ڶ��������ڵĵ�Ԫ������
        ISing=1
        Xie2(1) = Xi0
        Xie2(2) = Xie(2)
        Call SingElem2(Ge,He,Ind,Xie2,p2,knots2,KN2,weights2,PN2,cpts2,SXiN,ISing,tol)
        !
        ISing=2
        Xie2(1) = Xie(1)
        Xie2(2) = Xi0
        Call SingElem2(Ge2,He2,Ind2,Xie2,p2,knots2,KN2,weights2,PN2,cpts2,SXiN,ISing,tol)
        !
        if (any(ind /= ind2) /= 0) then
            print*, "�ڶ��������ڵĵ�Ԫ���������������, GeHe2.f90"
            read(*,*)
            stop
        end if
        Ge = Ge + Ge2
        He = He + He2
    elseif (l( 8) .or. l( 9)) then ! ����ڵ�һ����������ڣ���Ԫ�ڵڶ��������ڵķ�����
        Call NonSingElem2(Ge,He,ind,Np,C0,Xie,p2,knots2,KN2,weights2,PN2,cpts2)
    end if
    return
    end subroutine GeHe2