    subroutine collpoint(p,U,KN,W,WN,cpts,i,Xi0,RN,IEN,C0)
    implicit none
    ! 传递进来的变量声明
    integer(kind=4)::p,KN,WN,i
    integer(kind=4)::IEN(0:p,1:1)
    real(kind=8)::Xi0
    real(kind=8)::U(0:KN),W(0:WN),cpts(0:WN,1:2),RN(0:p,1:1),C0(1:2)
    ! 其他变量声明
    integer(kind=4)::j
    real(kind=8)::Xi(1:1)
    
    Xi=0.0D0
    do j=i+1,i+p,1
        Xi(1) = Xi(1) + U(j)
    end do
    Xi(1) = Xi(1) / p
    Call nrbasis(U,KN,Xi,1,W,WN,p,RN,IEN)
    Call RXicpts(RN,p,1,IEN,cpts,WN,2,C0)
    Xi0 = Xi(1) 
    
    return
    end subroutine collpoint