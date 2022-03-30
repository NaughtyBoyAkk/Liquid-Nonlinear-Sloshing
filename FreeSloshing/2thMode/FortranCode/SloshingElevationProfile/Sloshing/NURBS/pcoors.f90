    subroutine pcoors(p,U,KN,W,WN,cpts,Nc,Xi,XiN,Ci)
    implicit none
    ! 传递进来的变量声明
    integer(kind=4)::p,KN,WN,Nc,XiN
    real(kind=8)::U(0:KN),W(0:WN),Xi(1:XiN),cpts(0:WN,1:Nc),Ci(1:XiN,1:2)
    ! 其他进来的变量声明
    integer(kind=4)::IEN(0:p,1:XiN)
    real(kind=8)::RN(0:p,1:XiN)
    
    Call nrbasis(U,KN,Xi,XiN,W,WN,p,RN,IEN)
    Call RXicpts(RN,p,XiN,IEN,cpts,WN,Nc,Ci)

    return
    end subroutine pcoors