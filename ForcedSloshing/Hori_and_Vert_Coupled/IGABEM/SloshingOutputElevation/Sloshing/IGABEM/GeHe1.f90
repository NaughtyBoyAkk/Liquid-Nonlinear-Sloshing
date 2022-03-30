    subroutine GeHe1(Ge,He,Ind,Xi0,C0,Xie,p1,knots1,KN1,weights1,PN1,cpts1,l,Np,SXiN,tol)
    implicit none
    !���ݽ����ı�������
    integer(kind=4)::KN1,PN1,p1,Np,SXiN
    integer(kind=4)::Ind(0:p1)
    real(kind=8)::Xi0,tol
    real(kind=8)::knots1(0:KN1),cpts1(0:PN1,1:2),weights1(0:PN1)
    real(kind=8)::Ge(0:p1),He(0:p1),Xie(1:2),C0(1:2)
    logical::l(1:10)
    ! ������������
    integer(kind=4)::ISing
    integer(kind=4)::Ind2(0:p1)
    real(kind=8)::Xie2(1:2),Ge2(0:p1),He2(0:p1)
    !
    if (count(l) /= 1) then
        print*, "����������жϴ���, GeHe1.f90"
        read(*,*)
        stop
    elseif (l( 1)) then ! ��һ�������ڵĵ�һ�ڵ�����
        ISing=1
        Call SingElem2(Ge,He,Ind,Xie,p1,knots1,KN1,weights1,PN1,cpts1,SXiN,ISing,tol)   
    elseif (l( 3)) then ! ��һ�������ڵĵڶ��ڵ�����
        ISing=2
        Call SingElem2(Ge,He,Ind,Xie,p1,knots1,KN1,weights1,PN1,cpts1,SXiN,ISing,tol)
    elseif (l( 5)) then ! ��һ�������ڵĵ�Ԫ������
        ISing=1
        Xie2(1) = Xi0
        Xie2(2) = Xie(2)
        Call SingElem2(Ge,He,Ind,Xie2,p1,knots1,KN1,weights1,PN1,cpts1,SXiN,ISing,tol)
        !
        ISing=2
        Xie2(1) = Xie(1)
        Xie2(2) = Xi0
        Call SingElem2(Ge2,He2,Ind2,Xie2,p1,knots1,KN1,weights1,PN1,cpts1,SXiN,ISing,tol)
        if (any(ind /= ind2) /= 0) then
            print*, "��һ�������ڵĵ�Ԫ���������������, GeHe1.f90"
            read(*,*)
            stop
        end if
        !
        Ge = Ge + Ge2
        He = He + He2
    elseif (l( 7) .or. l(10)) then ! ����ڵ�һ����������ڣ���Ԫ�ڵ�һ�������ڵķ�����
        Call NonSingElem2(Ge,He,ind,Np,C0,Xie,p1,knots1,KN1,weights1,PN1,cpts1)
    end if
    return
    end subroutine GeHe1