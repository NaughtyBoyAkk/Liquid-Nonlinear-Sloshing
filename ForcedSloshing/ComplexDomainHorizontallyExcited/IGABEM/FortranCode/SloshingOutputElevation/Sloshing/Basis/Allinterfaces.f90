    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    interface
    !
    subroutine GeomPara(p,knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2)
    implicit none
    integer(kind=4)::p,KN1,PN1,KN2,PN2
    real(kind=8),allocatable::knots1(:),cpts1(:,:),weights1(:),knots2(:),cpts2(:,:),weights2(:)
    end subroutine GeomPara
    !
    subroutine UniqueKnots(uniknots,uniKN,knots,KN,tol)
    implicit none
    integer(kind=4)::KN,uniKN
    real(kind=8)::tol
    real(kind=8)::knots(0:KN) 
    real(kind=8),allocatable::uniknots(:)
    end subroutine UniqueKnots
    !
    
    end interface
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++