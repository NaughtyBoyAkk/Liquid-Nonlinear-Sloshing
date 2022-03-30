    subroutine Shell_nrbasis(U,KN,Xi,XiN,W,WN,p,RN,IEN)
    implicit none
    !
    integer(kind=4)::KN,XiN,WN,p
    integer(kind=4),intent(in)::IEN(0:p,1:XiN) ! 此处的IEN为指定要取出的索引，利用索引区分基函数在节点处的左右值
    real(kind=8)::U(0:KN),Xi(1:XiN),W(0:WN),RN(0:p,1:XiN)
    !
    integer(kind=4)::i
    Integer(kind=4)::IEN0(0:p,1:XiN),IEN0i(0:p), IENi(0:p)
    real(kind=8)::RNs(0:WN,1:XiN)
    !
    Call nrbasis(U,KN,Xi,XiN,W,WN,p,RN,IEN0)
    RNs = 0.0D0
    do i=1,XiN,1
        IEN0i = IEN0(0:p,i)
        RNs(IEN0i,i) = RN(0:p,i)
        IENi = IEN(0:p,i)
        RN(0:p,i) = RNs(IENi,i)
    end do
    return
    end subroutine Shell_nrbasis