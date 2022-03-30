    subroutine UniqueKnots(uniknots,uniKN,knots,KN,tol)
    implicit none
    integer(kind=4)::KN,uniKN
    real(kind=8)::tol
    real(kind=8)::knots(0:KN)
    real(kind=8),allocatable::uniknots(:)
    !������������
    integer(kind=4)::i,j
    real(kind=8)::knoti,knotj
    ! ��������
    j=0
    do i=0,KN-1,1
        knoti=knots(i)
        knotj=knots(i+1)
        if (abs(knoti-knotj) > tol) j=j+1
    end do
    ! �������飬����ֵ
    uniKN=j
    allocate(uniknots(0:j))
    uniknots(0)=knots(0)
    j=0
    do i=0,KN-1,1
        knoti=knots(i)
        knotj=knots(i+1)
        if (abs(knoti-knotj) > tol) then
            j=j+1
            uniknots(j)=knotj
        end if
    end do    
    return
    end subroutine UniqueKnots