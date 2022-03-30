    subroutine GetVN(VN,fileName)
    implicit none
    integer(kind=4)VN,ios1
    Character(len=1024)fileName
    character(len=1)Tempchar
    open(unit=666,file=fileName,form='formatted',access='sequential',status='old')
    VN=0
    do while(.True.)
        read(666,'(A1)',iostat=ios1)Tempchar
        if (ios1==0) then
            VN=VN+1
        else
            exit
        end if
    end do
    close(666)
    return
    end subroutine GetVN