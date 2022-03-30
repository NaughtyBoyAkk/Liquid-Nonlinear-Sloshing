    subroutine FindSpan(n,p,ui,U,KN,j)
    implicit none
    !传递进来的变量声明
    integer(kind=4)::n,p,KN,j
    real(kind=8)::ui
    real(kind=8)::U(0:KN)
    !其他变量声明
    integer(kind=4)::low,high,mid
    !**************************************
    if (ui < U(0) .or. ui > U(KN)) then
        print*, "u = ", ui, "should be located in U = [", U(0), ",", U(KN), "]"
        read(*,*)
        stop
    end if
    !
    if (ui==U(n+1)) then
        j=n
    else
        low=p
        high=n+1
        mid=INT((low+high)/2.0)
        do while (ui<U(mid) .OR. ui>=U(mid+1))
            if (ui<U(mid)) then
                high=mid
            else
                low=mid
            end if
            mid=INT((low+high)/2.0)
        end do
        j=mid
    end if
    return
    end subroutine FindSpan
