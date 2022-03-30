    subroutine DersBasisFuns(i,ui,p,dn,U,KN,ders)
    implicit none
    !传递进来的变量声明
    integer(kind=4)::i,p,dn,KN
    real(kind=8)::ui
    real(kind=8)::U(0:KN),ders(0:p,0:dn)
    !其他变量声明
    integer(kind=4)::j,r,k,s1,s2,j1,j2,rk,pk
    real(kind=8)::temp,saved,d
    real(kind=8)::ndu(0:p,0:p),a(0:p,0:1),left(0:KN),right(0:KN)
    !
    if (dn>p) then
        print*, "dn cannot be greater than p"
        read(*,*)
        stop
    end if
    ndu(0,0)=1.0D0
    do j=1,p,1
        left(j)=ui-U(i+1-j)
        right(j)=U(i+j)-ui
        saved=0.0D0
        do r=0,j-1,1
            ndu(j,r)=right(r+1)+left(j-r)
            temp=ndu(r,j-1)/ndu(j,r)
            ndu(r,j)=saved+right(r+1)*temp
            saved=left(j-r)*temp
        end do
        ndu(j,j)=saved
    end do
    do j=0,p,1
        ders(j,0)=ndu(j,p)
    end do
    !ders(0:p,0)=ndu(0:p,p)
    do r=0,p,1
        s1=0
        s2=1
        a(0,0)=1.0D0
        do k=1,dn,1
            d=0.0D0
            rk=r-k
            pk=p-k
            if (r>=k) then
                a(0,s2)=a(0,s1)/ndu(pk+1,rk)
                d=a(0,s2)*ndu(rk,pk)
            end if
            if (rk>=-1) then
                j1=1
            else
                j1=-rk
            end if
            if (r-1<=pk) then
                j2=k-1
            else
                j2=p-r
            end if
            do j=j1,j2,1
                a(j,s2)=(a(j,s1)-a(j-1,s1))/ndu(pk+1,rk+j)
                d=d+a(j,s2)*ndu(rk+j,pk)
            end do
            if (r<=pk) then
                a(k,s2)=-a(k-1,s1)/ndu(pk+1,r)
                d=d+a(k,s2)*ndu(r,pk)
            end if
            ders(r,k)=d
            j=s1
            s1=s2
            s2=j
        end do
    end do
    !
    r=p
    do k=1,dn,1
        do j=0,p,1
            ders(j,k)=ders(j,k)*dble(r)
        end do
        r=r*(p-k)
    end do
    return
    end subroutine DersBasisFuns