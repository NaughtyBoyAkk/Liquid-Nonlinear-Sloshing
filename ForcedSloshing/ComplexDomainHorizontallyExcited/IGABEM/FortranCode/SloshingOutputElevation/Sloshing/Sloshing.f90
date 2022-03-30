    program Sloshing
    implicit none
    include 'mkl_lapack.fi'
    include 'mkl_blas.fi'
    !
    real(kind=8),parameter::t0 = 0.0D0, t1 = 10.0D0, timestep = 0.012D0
    !
    real(kind=8),parameter::c1 = 0.5D0, c2 = 1.0D0 - c1, c3 = 0.5D0 / c2, ah = 2.0D-3, omegah = 5.5D0
    real(kind=8),parameter::tol = 1.0D-10, Accg = 9.81D0
    integer(kind=4),parameter::Np = 5, SXiN = 12, p = 3
    integer(kind=4)::KN1,KN2,PN1,PN2,EN1,EN2,p1,i
    integer(kind=4)::keypoints(1:4)
    real(kind=8)::AccX, ti
    real(kind=8),allocatable::knots1(:),cpts1(:,:),weights1(:),knots2(:),cpts2(:,:),weights2(:),nodes1(:),nodes2(:)
    !
    real(kind=8),allocatable::Cx(:,:),etafp(:),phifp(:),etacpts(:)
    real(kind=8),allocatable::dphidx(:),dphidz(:),detadx(:),dphidt(:),detadt(:)
    real(kind=8),allocatable::midetafp(:),midphifp(:)

    !
    include "./Basis/Allinterfaces.f90"
    ! 获取几何信息：节点向量，控制点坐标，控制点权重系数
    Call GeomPara(p,knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2)
    ! 获取非重复的节点向量
    Call UniqueKnots(nodes1,EN1,knots1,KN1,tol)
    Call UniqueKnots(nodes2,EN2,knots2,KN2,tol)
    !
    p1 = KN1 - PN1 - 1
    if (p1 /= p) then
        print*, "曲线阶次设置出错"
        read(*,*)
        stop
    end if
    ! keypoints 代表的是角点处控制点的编号，用来更新左右两侧的Wall的控制点
    keypoints(1) =  0
    keypoints(2) =  22
    keypoints(3) =  80
    keypoints(4) =  PN2
    ! 开始计算
    ! 初始条件
    allocate(Cx(0:EN1,1:2),etafp(0:EN1),phifp(0:EN1),etacpts(0:EN1),dphidx(0:EN1),dphidz(0:EN1),detadx(0:EN1))
    allocate(midetafp(0:EN1),midphifp(0:EN1))
    Call pcoors(p1,knots1,KN1,weights1,PN1,cpts1,2,nodes1,EN1+1,Cx)
    etafp = Cx(:,2)
    phifp = 0.0D0
    i = 0
    ti = t0
    open(666,file='.\results\elevation.dat',status='replace')
    do while (ti < t1)
        write(*,'(A,I10,A,F9.4)') " i=",i,"  ==>  ti= ",ti
        ! 输出波高
        write(666,'(F10.4,Es18.9)') ti, etafp(EN1)
        !
        Call PartialTerms(knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2,nodes1,EN1,nodes2,EN2, &
            & Np,SXIN,etafp,phifp,keypoints,dphidx,dphidz,detadx)
        AccX = -ah * omegah*omegah * sin( omegah * ti )
        detadt = dphidz - dphidx * detadx
        dphidt = dphidz * detadt - (dphidx * dphidx + dphidz * dphidz)/2.0D0 - etafp * Accg - Cx(:,1) *  AccX
        midetafp = etafp + c3 * timestep * detadt
        midphifp = phifp + c3 * timestep * dphidt
        Call PartialTerms(knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2,nodes1,EN1,nodes2,EN2, &
            & Np,SXIN,midetafp,midphifp,keypoints,dphidx,dphidz,detadx)
        AccX = -ah * omegah*omegah * sin( omegah * ( ti + c3 * timestep ) )
        detadt = dphidz - dphidx * detadx
        dphidt = dphidz * detadt - (dphidx * dphidx + dphidz * dphidz)/2.0D0 - midetafp * Accg - Cx(:,1) *  AccX
        etafp = (c3-c1)/c3 * etafp + c1/c3 * midetafp + c2 * timestep * detadt
        phifp = (c3-c1)/c3 * phifp + c1/c3 * midphifp + c2 * timestep * dphidt
        ! 进入下一时间
        ti = ti + timestep
        i=i+1
    end do
    close(666)
    print*," Completed! "
    !read(*,*)
    print *, 'Hello World'

    end program Sloshing

