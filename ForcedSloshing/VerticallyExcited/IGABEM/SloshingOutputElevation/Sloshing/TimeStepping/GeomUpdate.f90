    subroutine GeomUpdate(etacpts,cpts1,cpts2,PN1,PN2,keypoints)
    implicit none 
    !
    integer(kind=4)::PN1,PN2
    integer(kind=4)::keypoints(1:4)
    real(kind=8)::etacpts(0:PN1),cpts1(0:PN1,1:2),cpts2(0:PN2,1:2)
    !
    integer(kind=4)::k1,k2,k3,k4
    real(kind=8)::delta(1:2)
    real(kind=8),allocatable:: alpha1(:),deltay1(:),alpha2(:),deltay2(:)
    !
    k1 = keypoints(1)
    k2 = keypoints(2)
    k3 = keypoints(3)
    k4 = keypoints(4)
    !
    if (etacpts(PN1) < cpts2(k2,2)) then
        print*, "自由液面左侧位移过大"
        read(*,*)
        stop
    else if (etacpts(0) < cpts2(k3,2)) then 
        print*, "自由液面右侧位移过大"
        read(*,*)
        stop
    else
        continue
    end if 
    !
    delta(1) = etacpts(PN1) - cpts2(0,2)
    delta(2) = etacpts(0) - cpts2(PN2,2)
    
    allocate(alpha1(k1:k2),deltay1(k1:k2),alpha2(k3:k4),deltay2(k3:k4))
    alpha1(k1:k2) = (cpts2(k1:k2,2) - cpts2(k2,2)) / (cpts2(k1,2) - cpts2(k2,2))
    alpha2(k3:k4) = (cpts2(k3:k4,2) - cpts2(k3,2)) / (cpts2(k4,2) - cpts2(k3,2))
    
    deltay1 = alpha1 * delta(1) 
    deltay2 = alpha2 * delta(2) 
    
    cpts2(k1:k2,2) = cpts2(k1:k2,2) + deltay1(k1:k2)
    cpts2(k3:k4,2) = cpts2(k3:k4,2) + deltay2(k3:k4)
    cpts1(:,2) = etacpts
    
    end subroutine GeomUpdate