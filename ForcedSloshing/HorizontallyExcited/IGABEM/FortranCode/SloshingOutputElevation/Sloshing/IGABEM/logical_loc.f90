    subroutine logical_loc(Xi0,C0,Xie,Ce,point_in1,elem_in1,l,tol)
    implicit none
    ! ���ݽ����ı�������
    real(kind=8)::Xi0,tol
    real(kind=8)::Xie(1:2),C0(1:2),Ce(1:2,1:2)
    logical::point_in1,elem_in1
    logical::l(1:10)
    ! ���������ı�������
    real(kind=8)::Deltar(1:2)
    logical::s(1:7)
    !
    Deltar=sqrt((Ce(1:2,1)-C0(1))**2 + (Ce(1:2,2)-C0(2))**2)
    s( 1) = (Deltar(1) < tol) ! ����һ�ڵ�����
    s( 2) = (Deltar(2) < tol) ! ���ڶ��ڵ�����
    s( 3) = ((Xie(1)-Xi0) > tol) .or. ((Xi0-Xie(2)) > tol) ! ���λ�ڵ�Ԫ�ⲿ�����s[4]��s[5]ʹ��
    s( 4) = point_in1 .and. elem_in1 ! ���͵�Ԫ��λ�ڵ�һ��������
    s( 5) = (.not. point_in1) .and. (.not. elem_in1) ! ���͵�Ԫ��λ�ڵڶ���������
    s( 6) = point_in1 .and. (.not. elem_in1) ! ���λ�ڵ�һ�������ڣ���Ԫλ�ڵڶ���������
    s( 7) = (.not. point_in1) .and. elem_in1 ! ���λ�ڵڶ��������ڣ���Ԫλ�ڵ�һ��������
    l( 1) = s(1) .and. elem_in1 ! ��һ�������ڵĵ�һ�ڵ�����
    l( 2) = s(1) .and. (.not. elem_in1) ! �ڶ��������ڵĵ�һ�ڵ�����
    l( 3) = s(2) .and. elem_in1 ! ��һ�������ڵĵڶ��ڵ�����
    l( 4) = s(2) .and. (.not. elem_in1) ! �ڶ��������ڵĵڶ��ڵ�����
    l( 5) = (.not. s(1)) .and. (.not. s(2)) .and. (.not. s(3)) .and. s(4) ! ��һ�������ڵĵ�Ԫ������
    l( 6) = (.not. s(1)) .and. (.not. s(2)) .and. (.not. s(3)) .and. s(5) ! �ڶ��������ڵĵ�Ԫ������
    l( 7) = (.not. s(1)) .and. (.not. s(2)) .and. s(3) .and. s(4) ! ����ڵ�һ�������ڣ���Ԫ�ڵ�һ�������ڵķ�����
    l( 8) = (.not. s(1)) .and. (.not. s(2)) .and. s(3) .and. s(5) ! ����ڵڶ��������ڣ���Ԫ�ڵڶ��������ڵķ�����
    l( 9) = (.not. s(1)) .and. (.not. s(2)) .and. s(6) ! ����ڵ�һ�������ڣ���Ԫ�ڵڶ��������ڵķ�����
    l(10) = (.not. s(1)) .and. (.not. s(2)) .and. s(7) ! ����ڵڶ��������ڣ���Ԫ�ڵ�һ�������ڵķ�����

    return
    end subroutine logical_loc