    subroutine PartialTerms(knots1,cpts1,weights1,KN1,PN1,knots2,cpts2,weights2,KN2,PN2,nodes1,EN1,nodes2,EN2, &
        & Np,SXIN,etafp,phifp,keypoints,dphidx,dphidz,detadx,outi)
    implicit none
    !
    integer(kind=4)::KN1,KN2,PN1,PN2,EN1,EN2,Np,SXIN,outi
    integer(kind=4)::keypoints(1:4)
    real(kind=8)::knots1(0:KN1),knots2(0:KN2),cpts1(0:PN1,1:2),cpts2(0:PN2,1:2),weights1(0:PN1),weights2(0:PN2)
    real(kind=8)::nodes1(0:EN1),nodes2(0:EN2),etafp(0:EN1),phifp(0:EN1),dphidx(0:EN1),dphidz(0:EN1),detadx(0:EN1)
    !
    integer(kind=4)::p1,PNN,j
    integer(kind=4),allocatable::IENXi(:,:)
    real(kind=8),parameter::tol=1.0D-10
    real(kind=8),allocatable::phix(:),VeloNor(:),etacpts(:,:),phicpts(:,:),RNg(:,:),dRgdXI(:,:),dCgdXi(:,:),JXI(:),outn(:,:)
    real(kind=8),allocatable::dphidXi(:),dphidn(:)
    character(len=10)::timetemp
    character(len=1024)::dataf
    !
    p1=KN1-PN1-1
    allocate(etacpts(0:PN1,1:1),phicpts(0:PN1,1:1))
    ! ����Һ�沨�ߵĿ��Ƶ�
    Call CompuCpts(p1,knots1,Weights1,cpts1,KN1,PN1,Nodes1,EN1,Np,etafp,etacpts)
    ! ���¼�����״
    Call GeomUpdate(etacpts,cpts1,cpts2,PN1,PN2,keypoints)
    ! ! �߽�����-��������Һ����phi�߽������Ŀ��Ƶ�
    Call CompuCpts(p1,knots1,Weights1,cpts1,KN1,PN1,Nodes1,EN1,Np,phifp,phicpts)
    !***********************************************************************************
    ! Laplace 2D ���
    PNN=PN1+PN2+1 ! NN������ܿ��Ƶ���Ŀ�Ǵ�0��ʼ�Ƶģ�����������PN1+PN2+1��������PN1+PN2+2
    allocate(phix(0:PNN),VeloNor(0:PNN))
    phix(0:PN1) = phicpts(0:PN1,1)
    phix(PN1+1:PNN) = 0.0D0
    VeloNor = 0.0D0
    Call Laplace2D(knots1,knots2,Weights1,Weights2,cpts1,cpts2,nodes1,nodes2, & 
        & KN1,KN2,PN1,PN2,EN1,EN2,Np,SxiN,phix,VeloNor,PNN)
    !����ٶ��ƿ��Ƶ�
    if (Outi> -1) then
        write(timetemp,'(I8)') Outi
        dataf = '.\results\phicpts'//trim(adjustl(timetemp))//'.phi'
        open(666,file=dataf,status='replace')
        do j=0,PNN,1
            write(666,'(2Es18.9)') phix(j)
        end do
        close(666)
        !
        dataf = '.\results\vcpts'//trim(adjustl(timetemp))//'.vel'
        open(666,file=dataf,status='replace')
        do j=0,PNN,1
            write(666,'(2Es18.9)') VeloNor(j)
        end do
        close(666)
        !
    end if
    ! ����ƫ����
    allocate(outn(0:EN1,1:2),dphidXi(0:EN1),dphidn(0:EN1))
    allocate(RNg(0:p1,0:EN1),dRgdXI(0:p1,0:EN1),IENXi(0:p1,0:EN1),dCgdXi(0:EN1,1:2),JXI(0:EN1))
    Call Ders1nrbasis(knots1,KN1,nodes1,EN1+1,Weights1,PN1,p1,dRgdXI,IENXI)
    Call RXicpts(dRgdXI,p1,EN1+1,IENXi,cpts1,PN1,2,dCgdXi)
    JXi=sqrt(dCgdXi(0:EN1,1)*dCgdXi(0:EN1,1)+dCgdXi(0:EN1,2)*dCgdXi(0:EN1,2))
    Call nrbasis(knots1,KN1,nodes1,EN1+1,Weights1,PN1,p1,RNg,IENXI)
    outn(0:EN1,1) =  dCgdXi(0:EN1,2) / JXi
    outn(0:EN1,2) = -dCgdXi(0:EN1,1) / JXi
    JXi = dCgdXi(0:EN1,1) * outn(0:EN1,2) - dCgdXi(0:EN1,2) * outn(0:EN1,1) !JXi ��������
    Call RXicpts(dRgdXI,p1,EN1+1,IENXi,phix(0:PN1),PN1,1,dphidXi)
    Call RXicpts(RNg,p1,EN1+1,IENXi,VeloNor(0:PN1),PN1,1,dphidn)
    dphidx = ( outn(0:EN1,2) * dphidXi - dCgdXi(0:EN1,2) * dphidn ) / JXi
    dphidz = ( dCgdXi(0:EN1,1) * dphidn - outn(0:EN1,1) * dphidXI ) / JXi
    detadx = dCgdXi(0:EN1,2) / dCgdXi(0:EN1,1)
    !
    dphidx(0) = 0.0D0 
    dphidx(EN1) = 0.0D0 
    return
    end subroutine PartialTerms
