    subroutine logical_loc(Xi0,C0,Xie,Ce,point_in1,elem_in1,l,tol)
    implicit none
    ! 传递进来的变量声明
    real(kind=8)::Xi0,tol
    real(kind=8)::Xie(1:2),C0(1:2),Ce(1:2,1:2)
    logical::point_in1,elem_in1
    logical::l(1:10)
    ! 其他进来的变量声明
    real(kind=8)::Deltar(1:2)
    logical::s(1:7)
    !
    Deltar=sqrt((Ce(1:2,1)-C0(1))**2 + (Ce(1:2,2)-C0(2))**2)
    s( 1) = (Deltar(1) < tol) ! 配点第一节点奇异
    s( 2) = (Deltar(2) < tol) ! 配点第二节点奇异
    s( 3) = ((Xie(1)-Xi0) > tol) .or. ((Xi0-Xie(2)) > tol) ! 配点位于单元外部，配合s[4]和s[5]使用
    s( 4) = point_in1 .and. elem_in1 ! 配点和单元均位于第一条曲线内
    s( 5) = (.not. point_in1) .and. (.not. elem_in1) ! 配点和单元均位于第二条曲线内
    s( 6) = point_in1 .and. (.not. elem_in1) ! 配点位于第一条曲线内，单元位于第二条曲线内
    s( 7) = (.not. point_in1) .and. elem_in1 ! 配点位于第二条曲线内，单元位于第一条曲线内
    l( 1) = s(1) .and. elem_in1 ! 第一条曲线内的第一节点奇异
    l( 2) = s(1) .and. (.not. elem_in1) ! 第二条曲线内的第一节点奇异
    l( 3) = s(2) .and. elem_in1 ! 第一条曲线内的第二节点奇异
    l( 4) = s(2) .and. (.not. elem_in1) ! 第二条曲线内的第二节点奇异
    l( 5) = (.not. s(1)) .and. (.not. s(2)) .and. (.not. s(3)) .and. s(4) ! 第一条曲线内的单元内奇异
    l( 6) = (.not. s(1)) .and. (.not. s(2)) .and. (.not. s(3)) .and. s(5) ! 第二条曲线内的单元内奇异
    l( 7) = (.not. s(1)) .and. (.not. s(2)) .and. s(3) .and. s(4) ! 配点在第一条曲线内，单元在第一条曲线内的非奇异
    l( 8) = (.not. s(1)) .and. (.not. s(2)) .and. s(3) .and. s(5) ! 配点在第二条曲线内，单元在第二条曲线内的非奇异
    l( 9) = (.not. s(1)) .and. (.not. s(2)) .and. s(6) ! 配点在第一条曲线内，单元在第二条曲线内的非奇异
    l(10) = (.not. s(1)) .and. (.not. s(2)) .and. s(7) ! 配点在第二条曲线内，单元在第一条曲线内的非奇异

    return
    end subroutine logical_loc