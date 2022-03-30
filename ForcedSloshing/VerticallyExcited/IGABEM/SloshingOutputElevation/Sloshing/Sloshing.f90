    program Sloshing
    implicit none
    include 'mkl_lapack.fi'
    include 'mkl_blas.fi'
    !
    real(kind=8),parameter::t0 = 0.0D0, t1 = 20.0D0, timestep = 0.005D0
    !
    real(kind=8),parameter::c1 = 0.5D0, c2 = 1.0D0 - c1, c3 = 0.5D0 / c2
    real(kind=8),parameter::tol = 1.0D-10, Accg = 9.81D0, tankL = 1.0D0, tankH = 0.5D0
    integer(kind=4),parameter::Np = 5, SXiN = 10, p = 3
    integer(kind=4)::KN1,KN2,PN1,PN2,EN1,EN2,p1
    integer(kind=4)::keypoints(1:4)
    real(kind=8)::AccZ, ti, Amp, waveNum, omegaf, Omega1, omegav, kv, av, epsilon, wavelen, analytical_eta
    real(kind=8),allocatable::knots1(:),cpts1(:,:),weights1(:),knots2(:),cpts2(:,:),weights2(:),nodes1(:),nodes2(:)
    !
    real(kind=8),allocatable::Cx(:,:),etafp(:),phifp(:),etacpts(:)
    real(kind=8),allocatable::dphidx(:),dphidz(:),detadx(:),dphidt(:),detadt(:)
    real(kind=8),allocatable::midetafp(:),midphifp(:)
    !
    include "./Basis/Allinterfaces.f90"
    !
    ! ������Ĳ���
    wavelen = 2.0D0 * tankL
    waveNum = 2.0D0 * acos(-1.0D0) / wavelen
    omegaf = sqrt( Accg * waveNum * tanh(waveNum * tankH) )
    Omega1 = 1.253D0
    omegav = omegaf / Omega1
    kv = 0.5D0
    av = Accg * kv / omegav**2 
    !    
    epsilon = 0.0014D0
    Amp = Accg * epsilon / omegaf**2
    ! ��ȡ������Ϣ���ڵ����������Ƶ����꣬���Ƶ�Ȩ��ϵ��
    Call GeomPara(p,knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2)
    ! ��ȡ���ظ��Ľڵ�����
    Call UniqueKnots(nodes1,EN1,knots1,KN1,tol)
    Call UniqueKnots(nodes2,EN2,knots2,KN2,tol)
    !
    p1 = KN1 - PN1 - 1
    if (p1 /= p) then
        print*, "���߽״����ó���"
        read(*,*)
        stop
    end if
    ! keypoints ������ǵڶ����ߣ������ڣ��Ľǵ��λ�ã������������������Wall�Ŀ��Ƶ�
    keypoints(1) =  0
    keypoints(2) =  17
    keypoints(3) =  49
    keypoints(4) =  PN2
    ! ��ʼ����
    ! ��ʼ����
    allocate(Cx(0:EN1,1:2),etafp(0:EN1),phifp(0:EN1),etacpts(0:EN1),dphidx(0:EN1),dphidz(0:EN1),detadx(0:EN1))
    allocate(midetafp(0:EN1),midphifp(0:EN1))
    Call pcoors(p1,knots1,KN1,weights1,PN1,cpts1,2,nodes1,EN1+1,Cx)
    etafp = Cx(:,2)  ! ����Һ���ϳ�ʼ���ߵĿ��Ƶ�
    phifp = 0.0D0    ! ����Һ���ϳ�ʼ�ٶ��ƵĿ��Ƶ�
    !
    ti = t0
    open(666,file='.\results\data.dat',status='replace')
    do while (ti < t1)
        write(*,'(A,F9.4)') " ti= ",ti
        ! �������
        write(666,'(2Es18.9)') ti*omegaf, etafp(EN1)/Amp !���ǵ�
        !
        Call PartialTerms(knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2,nodes1,EN1,nodes2,EN2, &
            & Np,SXIN,etafp,phifp,keypoints,dphidx,dphidz,detadx)
        AccZ = -av * omegav*omegav * cos( omegav * ti )
        detadt = dphidz - dphidx * detadx
        dphidt = dphidz * detadt - (dphidx * dphidx + dphidz * dphidz)/2.0D0 - etafp * (Accg + AccZ) 
        midetafp = etafp + c3 * timestep * detadt
        midphifp = phifp + c3 * timestep * dphidt
        Call PartialTerms(knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2,nodes1,EN1,nodes2,EN2, &
            & Np,SXIN,midetafp,midphifp,keypoints,dphidx,dphidz,detadx)
        AccZ = -av * omegav*omegav * cos( omegav * ( ti + c3 * timestep ) )
        detadt = dphidz - dphidx * detadx
        dphidt = dphidz * detadt - (dphidx * dphidx + dphidz * dphidz)/2.0D0 - midetafp * (Accg + AccZ) 
        etafp = (c3-c1)/c3 * etafp + c1/c3 * midetafp + c2 * timestep * detadt
        phifp = (c3-c1)/c3 * phifp + c1/c3 * midphifp + c2 * timestep * dphidt
        ti = ti + timestep
    end do
    close(666)
    print*," Completed! "
    !read(*,*)
    print *, 'Hello World'

    end program Sloshing

