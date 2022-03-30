    program Sloshing
    implicit none
    include 'mkl_lapack.fi'
    include 'mkl_blas.fi'
    !
    real(kind=8),parameter::t0 = 0.0D0, t1 = 4.0D0, timestep = 0.012D0
    !
    real(kind=8),parameter::c1 = 0.5D0, c2 = 1.0D0 - c1, c3 = 0.5D0 / c2, ah = 0.0D0, omegah = 5.5D0
    real(kind=8),parameter::tol = 1.0D-10, Accg = 9.81D0, tankL = 0.9D0, tankH = 0.6D0
    integer(kind=4),parameter::Np = 5, SXiN = 10, p = 3
    integer(kind=4)::KN1,KN2,PN1,PN2,EN1,EN2,p1
    integer(kind=4)::keypoints(1:4)
    real(kind=8)::AccX, ti, Amp, waveNum, omegaf, wavelen, analytical_eta
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
    ! keypoints 代表的是角点的位置，用来更新左右两侧的Wall的控制点
    keypoints(1) =  0
    keypoints(2) =  22
    keypoints(3) =  54
    keypoints(4) =  PN2
    ! 开始计算
    ! 初始条件
    allocate(Cx(0:EN1,1:2),etafp(0:EN1),phifp(0:EN1),etacpts(0:EN1),dphidx(0:EN1),dphidz(0:EN1),detadx(0:EN1))
    allocate(midetafp(0:EN1),midphifp(0:EN1))
    Call pcoors(p1,knots1,KN1,weights1,PN1,cpts1,2,nodes1,EN1+1,Cx)
    etafp = Cx(:,2)  ! 自由液面上初始波高的控制点
    phifp = 0.0D0    ! 自由液面上初始速度势的控制点
    ! 解析解的参数
    wavelen = 2.0D0 * tankL
    Amp = Cx(EN1,2)   !左侧角点
    waveNum = 2.0D0 * acos(-1.0D0) / wavelen
    omegaf = sqrt( Accg * waveNum * tanh(waveNum * tankH) )
    !
    ti = t0
    open(666,file='.\results\elevation.dat',status='replace')
    do while (ti < t1)
        write(*,'(A,F9.4)') " ti= ",ti
        analytical_eta = Amp * cos(waveNum*(Cx(EN1,1)+tankL/2.0D0)) * cos(omegaf*ti) !左侧角点
        ! 输出波高
        write(666,'(F10.4,2Es18.9)') ti, analytical_eta, etafp(EN1) !左侧角点
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
        ti = ti + timestep
    end do
    close(666)
    print*," Completed! "
    !read(*,*)
    print *, 'Hello World'

    end program Sloshing

