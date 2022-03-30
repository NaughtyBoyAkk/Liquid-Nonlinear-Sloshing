    subroutine BasisFuns(i,ui,p,U,KN,BN)
    implicit none
    !传递进来的变量声明
    integer(kind=4)::i,p,KN
    real(kind=8)::ui
    real(kind=8)::U(0:KN),BN(0:p)
    !其他变量声明
    integer(kind=4)::j,r
    real(kind=8)::saved,temp
    real(kind=8)::left(0:KN),right(0:KN)
    !
    BN(0)=1.0D0
    do j=1,p,1
        left(j)=ui-U(i+1-j)
        right(j)=U(i+j)-ui
        saved=0.0D0
        do r=0,j-1,1
            temp=BN(r)/(right(r+1)+left(j-r))
            BN(r)=saved+right(r+1)*temp
            saved=left(j-r)*temp
        end do
        BN(j)=saved
    end do
    return
    end subroutine BasisFuns