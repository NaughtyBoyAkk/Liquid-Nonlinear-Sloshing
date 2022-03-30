    subroutine nrbasis(U,KN,Xi,XiN,W,WN,p,RN,IEN)
    implicit none
    !传递进来的变量声明
    integer(kind=4)::KN,XiN,WN,p
    integer(kind=4)::IEN(0:p,1:XiN)
    real(kind=8)::U(0:KN),Xi(1:XiN),W(0:WN),RN(0:p,1:XiN)
    !其他变量声明
    integer(kind=4)::i,j,k
    integer(kind=4)::IENi(0:p)
    real(kind=8)::BN(0:p,1:XiN),NTW(1:XiN)
    !
    Call bsbasis(U,KN,p,Xi,XiN,BN,IEN)
    NTW=0.0D0
    do i=1,XiN,1
        do j=0,p,1
            k=IEN(j,i)
            NTW(i)=NTW(i)+BN(j,i)*W(k)
        end do
        IENi=IEN(0:p,i)
        RN(0:p,i)=BN(0:p,i)*W(IENi)/NTW(i)
    end do
    return
    end subroutine nrbasis