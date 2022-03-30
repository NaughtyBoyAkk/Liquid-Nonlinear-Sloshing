    subroutine Laplace2D(knots1,knots2,Weights1,Weights2,cpts1,cpts2,nodes1,nodes2, & 
        & KN1,KN2,PN1,PN2,EN1,EN2,Np,SXiN,phix,Veloy,NN)
    implicit none
    !
    integer(kind=4)::KN1,KN2,PN1,PN2,EN1,EN2,Np,SxIn,NN
    real(kind=8)::knots1(0:KN1),knots2(0:KN2),Weights1(0:PN1),Weights2(0:PN2),cpts1(0:PN1,1:2),cpts2(0:PN2,1:2)
    real(kind=8)::nodes1(0:EN1),nodes2(0:EN2),phix(1:NN+1),Veloy(1:NN+1)
    !
    integer(kind=4)::p1,p2,EN,i,j,k
    integer(kind=4),allocatable::Ind(:),IEN0(:,:),ipiv(:)
    real(kind=8),parameter::tol=1.0D-10
    real(kind=8)::Xi0,Calpha
    real(kind=8)::C0(1:2),Xie(1:2),Ce(1:2,1:2)
    real(kind=8),allocatable::RN0(:,:)
    real(kind=8),allocatable::Ge(:),He(:),G(:,:),H(:,:),H1temp(:,:),sol(:)
    logical::point_in1,elem_in1
    logical::l(1:10)
    !
    p1=KN1-PN1-1
    p2=KN2-PN2-1
    ! 开始计算系数矩阵
    !
    EN=EN1+EN2
    allocate(RN0(0:p1,1:1),IEN0(0:p1,1:1),Ge(0:p2),He(0:p1),Ind(0:p1),G(1:NN+1,1:NN+1),H(1:NN+1,1:NN+1))
    G = 0.0D0
    H = 0.0D0
    do i=0,NN,1    ! 节点循环
        point_in1 = (i <= PN1)
        if (point_in1) then
            Call collpoint(p1,knots1,KN1,weights1,PN1,cpts1,i,Xi0,RN0,IEN0,C0)
            ind = IEN0(0:p1,1) + 1
        else
            k=i-PN1-1
            Call collpoint(p2,knots2,KN2,weights2,PN2,cpts2,k,Xi0,RN0,IEN0,C0)
            ind = IEN0(0:p1,1) + 1 + PN1 + 1
        end if
        ! Jump terms
        Call Jumpterms(Calpha,i,XI0,p1,knots1,Weights1,cpts1,p2,knots2,Weights2,cpts2,KN1,KN2,PN1,PN2,tol)
        H(i+1,ind) = H(i+1,ind) + Calpha * RN0(0:p1,1)
        ! 单元积分
        do j=1,EN,1    ! 单元循环
            elem_in1 = (j <= EN1)
            if (elem_in1) then
                Xie(1) = nodes1(j-1)
                Xie(2) = nodes1(j)
                Call pcoors(p1,knots1,KN1,weights1,PN1,cpts1,2,Xie,2,Ce)
                ! 判断节点和单元位置关系
                Call logical_loc(Xi0,C0,Xie,Ce,point_in1,elem_in1,l,tol)
                Call GeHe1(Ge,He,Ind,Xi0,C0,Xie,p1,knots1,KN1,weights1,PN1,cpts1,l,Np,SXiN,tol)
            else
                k=j-EN1
                Xie(1) = nodes2(k-1)
                Xie(2) = nodes2(k)
                Call pcoors(p2,knots2,KN2,weights2,PN2,cpts2,2,Xie,2,Ce)
                ! 判断节点和单元位置关系
                Call logical_loc(Xi0,C0,Xie,Ce,point_in1,elem_in1,l,tol)
                Call GeHe2(Ge,He,Ind,Xi0,C0,Xie,p2,knots2,KN2,weights2,PN2,cpts2,l,Np,SXiN,tol)
                Ind = Ind + PN1 + 1
            end if
            G(i+1,Ind) = G(i+1,Ind) + Ge
            H(i+1,Ind) = H(i+1,Ind) + He
        end do
    end do
    ! Boundary Conditions
    allocate(H1temp(1:NN+1,1:PN1+1),sol(1:NN+1),ipiv(1:NN+1))
    !
    H1temp = H(1:NN+1,1:PN1+1)
    H(1:NN+1,1:PN1+1) = -G(1:NN+1,1:PN1+1)
    G(1:NN+1,1:PN1+1) = -H1temp
    Call dgemv('N', NN+1, NN+1, 1.0D0, G, NN+1, phix, 1, 0.0D0, sol, 1)
    ! 求解
    Call dgetrf(NN+1,NN+1,H,NN+1,ipiv,i)
    if (i /= 0) then
        print*, "LU分解出错, Sloshing.f90"
        read(*,*)
        stop
    end if
    Call dgetrs('N',NN+1,1,H,NN+1,ipiv,sol,NN+1,j)
    if (j /= 0) then
        print*, "求解出错, Sloshing.f90"
        read(*,*)
        stop
    end if
    Veloy(1:PN1+1) = sol(1:PN1+1)
    Veloy(PN1+2:NN+1) = phix(PN1+2:NN+1)
    phix(PN1+2:NN+1) = sol(PN1+2:NN+1)
    return
    end subroutine Laplace2D