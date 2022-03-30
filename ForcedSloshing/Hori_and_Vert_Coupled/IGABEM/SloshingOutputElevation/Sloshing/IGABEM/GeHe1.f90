    subroutine GeHe1(Ge,He,Ind,Xi0,C0,Xie,p1,knots1,KN1,weights1,PN1,cpts1,l,Np,SXiN,tol)
    implicit none
    !传递进来的变量声明
    integer(kind=4)::KN1,PN1,p1,Np,SXiN
    integer(kind=4)::Ind(0:p1)
    real(kind=8)::Xi0,tol
    real(kind=8)::knots1(0:KN1),cpts1(0:PN1,1:2),weights1(0:PN1)
    real(kind=8)::Ge(0:p1),He(0:p1),Xie(1:2),C0(1:2)
    logical::l(1:10)
    ! 其他变量声明
    integer(kind=4)::ISing
    integer(kind=4)::Ind2(0:p1)
    real(kind=8)::Xie2(1:2),Ge2(0:p1),He2(0:p1)
    !
    if (count(l) /= 1) then
        print*, "配点奇异性判断错误, GeHe1.f90"
        read(*,*)
        stop
    elseif (l( 1)) then ! 第一条曲线内的第一节点奇异
        ISing=1
        Call SingElem2(Ge,He,Ind,Xie,p1,knots1,KN1,weights1,PN1,cpts1,SXiN,ISing,tol)   
    elseif (l( 3)) then ! 第一条曲线内的第二节点奇异
        ISing=2
        Call SingElem2(Ge,He,Ind,Xie,p1,knots1,KN1,weights1,PN1,cpts1,SXiN,ISing,tol)
    elseif (l( 5)) then ! 第一条曲线内的单元内奇异
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
            print*, "第一条曲线内的单元内奇异计算有问题, GeHe1.f90"
            read(*,*)
            stop
        end if
        !
        Ge = Ge + Ge2
        He = He + He2
    elseif (l( 7) .or. l(10)) then ! 配点在第一或二条曲线内，单元在第一条曲线内的非奇异
        Call NonSingElem2(Ge,He,ind,Np,C0,Xie,p1,knots1,KN1,weights1,PN1,cpts1)
    end if
    return
    end subroutine GeHe1