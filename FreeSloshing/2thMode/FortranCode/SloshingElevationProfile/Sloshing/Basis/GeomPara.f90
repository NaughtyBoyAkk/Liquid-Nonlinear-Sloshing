    subroutine GeomPara(p,knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2)
    implicit none
    integer(kind=4)::p,KN1,PN1,KN2,PN2
    real(kind=8),allocatable::knots1(:),cpts1(:,:),weights1(:),knots2(:),cpts2(:,:),weights2(:)
    !其他变量声明
    integer(kind=4)::i
    real(kind=8)::temp
    character(len=6)::pc
    character(len=1024)::fileknot1,filecpt1,fileknot2,filecpt2
    !
    write(pc,'(I4)') p
    fileknot1 =  "./Geom/knots1_p"//trim(adjustl(pc))//".txt"
    filecpt1  =  "./Geom/cpts1_p"//trim(adjustl(pc))//".txt"
    fileknot2 =  "./Geom/knots2_p"//trim(adjustl(pc))//".txt"
    filecpt2  =  "./Geom/cpts2_p"//trim(adjustl(pc))//".txt"
    Call GetVN(KN1,fileknot1)
    Call GetVN(PN1,filecpt1)
    Call GetVN(KN2,fileknot2)
    Call GetVN(PN2,filecpt2)
    KN1=KN1+1
    PN1=PN1-1
    KN2=KN2+1
    PN2=PN2-1
    allocate(knots1(0:KN1),cpts1(0:PN1,1:2),weights1(0:PN1),knots2(0:KN2),cpts2(0:PN2,1:2),weights2(0:PN2))
    !几何1
    open(unit=666,file=fileknot1,form='formatted',access='sequential',status='old')
    do i=1,KN1-1,1
        read(666,*)knots1(i)
    end do
    close(666)
    knots1(0)=knots1(1)
    knots1(KN1)=knots1(KN1-1)
    !
    open(unit=777,file=filecpt1,form='formatted',access='sequential',status='old')
    do i=0,PN1,1
        read(777,*)cpts1(i,1),cpts1(i,2),temp,weights1(i)
    end do
    close(777)
    !几何2
    open(unit=666,file=fileknot2,form='formatted',access='sequential',status='old')
    do i=1,KN2-1,1
        read(666,*)knots2(i)
    end do
    close(666)
    knots2(0)=knots2(1)
    knots2(KN2)=knots2(KN2-1)
    !
    open(unit=777,file=filecpt2,form='formatted',access='sequential',status='old')
    do i=0,PN2,1
        read(777,*)cpts2(i,1),cpts2(i,2),temp,weights2(i)
    end do
    close(777)

    return
    end subroutine GeomPara