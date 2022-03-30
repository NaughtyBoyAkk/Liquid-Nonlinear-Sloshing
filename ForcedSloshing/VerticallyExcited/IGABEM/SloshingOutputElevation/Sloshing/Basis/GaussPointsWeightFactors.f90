    subroutine GaussPointsWeightFactors(x1,x2,x,w,n)
    implicit none
    real(kind=8) x1,x2,EPS,pi,p1,p2,p3,pp,xl,xm,z,z1
    real(kind=8)x(n),w(n)
    integer(kind=4) i,j,m,n
    eps=1.0d-14
    !对于双精度(kind=8)存储的实型来说，这个精度已经很高，如果想把精度提高到1.0d-16以上
    !那么得将实型的kind值再提高，令其大于8
    pi=acos(-1.0d0)
    !EPS is the relative precision.
    !Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
    !arrays x(1:n) and w(1:n) of length n, containing the abscissas and weights of the Gauss-
    !Legendre n-point quadrature formula.
    !High precision is a good idea for this routine.
    m=(n+1)/2 !The roots are symmetric in the interval, so we
    xm=0.5d0*(x2+x1) !only have to nd half of them.
    xl=0.5d0*(x2-x1)
    do i=1,m,1 !Loop over the desired roots.
        z=cos(pi*(dble(i)-0.25d0)/(dble(n)+0.5d0))
        !Starting with the above approximation to the ith root, we enter the main loop of re-
        !nement by Newton's method.
1       continue
        p1=1.0d0
        p2=0.0d0
        do j=1,n,1 !Loop up the recurrence relation to get the Leg
            p3=p2   !endre polynomial evaluated at z.
            p2=p1
            p1=((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/dble(j)
        end do
        !p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by
        !a standard relation involving also p2, the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.0d0)
        z1=z
        z=z1-p1/pp !Newton's method.
        if (abs(z-z1)>EPS) goto 1
        x(i)=xm-xl*z !Scale the root to the desired interval,
        x(n+1-i)=xm+xl*z !and put in its symmetric counterpart.
        w(i)=2.0d0*xl/((1.d0-z*z)*pp*pp) !Compute the weight
        w(n+1-i)=w(i) !and its symmetric counterpart.
    end do
    return
    end