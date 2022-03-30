    subroutine bsbasis(U,KN,p,Xi,XiN,BN,IEN)
    implicit none
    !传递进来的变量声明
    integer(kind=4)::KN,p,XiN
    integer(kind=4)::IEN(0:p,1:XiN)
    real(kind=8)::U(0:KN),Xi(1:XiN),BN(0:p,1:XiN)
    !其他变量声明
    integer(kind=4)::i,j,n,k
    real(kind=8)::xii
    real(kind=8)::BNi(0:p)
    !
    n=KN-p-1
    do i=1,XiN,1
        xii=Xi(i)
        Call FindSpan(n,p,xii,U,KN,j)
        Call BasisFuns(j,xii,p,U,KN,BNi)
        BN(0:p,i)=BNi
        do k=0,p,1
            IEN(p-k,i)=j-k
        end do
    end do
    return
    end subroutine bsbasis