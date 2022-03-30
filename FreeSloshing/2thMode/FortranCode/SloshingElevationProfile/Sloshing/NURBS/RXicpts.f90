    subroutine RXicpts(RN,p,XiN,IEN,cpts,PN,Nc,Val)
    implicit none
    !传递进来的参数声明
    integer(kind=4)::p,XiN,PN,Nc
    integer(kind=4)::IEN(0:p,1:XiN)
    real(kind=8)::RN(0:p,1:XiN),cpts(0:PN,1:Nc),Val(1:XiN,1:Nc)
    !其他变量声明
    integer(kind=4)::i,j,k,Nci
    !
    Val = 0.0D0
    do Nci=1,Nc,1
        do i=1,XiN,1
            do j=0,p,1
                k=IEN(j,i)
                Val(i,Nci)=Val(i,Nci)+RN(j,i)*cpts(k,Nci)
            end do
        end do
    end do
    return
    end subroutine RXicpts
