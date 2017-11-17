C----------------------------------------------------------------------
C  This code produces a grid in x and Q with the values of the PDFs.
C  The code is based on the tabulation_example code of HOPPET
C 
C  Needs to be linked to cernlib
C----------------------------------------------------------------------
      program xpdfs_grid
      implicit none
      !--- variables defining the grid and evolution parameters
      double precision dy, ymax, Qmin, Qmax, dlnlnQ
      integer          nloop, order
      !---------------
      external         f1model
      double precision Q, ourpdf(-6:6)
      double precision asQ0, Q0, xmu ,Qmodel, asQmodel
      double precision hoppetalphas
      double precision xmin, xmax, xstep, xi
      double precision QQmin, QQmax, QQstep, QQi
      double precision Qin(1:3), Qfin(1:3)
      integer          nxstep, nQQstep
      integer          ix, iQQ, scheme, i, ithre
      integer          ipdf
      double precision cthre, bthre, tthre

      common/pdftype/ipdf

      common/alphas/Qmodel

      ! start the dglap evolution/convolution package 
      ymax  = 12d0      ! max value of ln 1/x
      dy    = 0.1d0     ! the internal grid spacing (smaller->higher accuarcy)
                        ! 0.1 should provide at least 10^{-3} accuracy 
      Qmin  = dsqrt(0.3d0)     ! smallest Q value in tabulation
      Qmax  = 10d0       ! largest Q value in tabulation
      dlnlnQ = dy/4     ! tabulation spacing in dlnlnQ (dy/4 recommended)
      nloop  = 1        ! the number of loops to initialise (max=3!)
      order  = -6       ! numerical interpolation order (-6 is a good choice)

      ! change the ipdf to select different pdfs
      ipdf = 2          ! 1=f1, 2=f1Tp1, 3=g1, 4=g1T1, 5=h1, 
                        ! 6=h1p1, 7=h1Lp1, 8=h1Tp1

c- the number of step is the following times 3, because the grid is done in 
c- 3 parts 
      nQQstep= 15

c- here fix the values of the charm, bottom, and top masses
*      cthre=1.414213563d0
      cthre=1.4d0
      bthre=4.5d0
      tthre=175d0

c- here select the values of the min and max z and number of steps
      xmin=0.01
      xmax=0.9
      nxstep= 50
      xstep=(xmax- xmin)/dble(nxstep-1)

      ! set the initial scale and the coupling there (in general the PDF
      ! and coupling may be specified at different scales -- this is not
      ! done here) and 
      Q0   = 91.187
      asQ0 = 0.125
*      Q0   = dsqrt(0.3d0)
*      asQ0 = 0.706862

      ! the ratio xmu = mu_F/mu_R to be used in the evolution.
      xmu  = 1.0d0

c- choose the initial Q of model or parametrization
      Qmodel=dsqrt(0.3d0)


      if(ipdf.eq.1.or.ipdf.eq.2) then
        scheme = 1        ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar, 
                          ! 4=frag, 5=tranPol-LO
      else if(ipdf.eq.3.or.ipdf.eq.4) then
        scheme = 3
      else
        scheme = 5     
      end if

      ! call this once at the beginning of your program
      call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,
     $                         scheme)
      write(6,'(a)') "Splitting functions initialised!"
      
      ! tell hoppet to use a variable flavour number scheme with
      ! the following c,b,t quark masses
      call hoppetSetVFN(cthre,bthre,tthre)


      ! carry out the evolution to create a tabulation, corresponding
      ! to the initial condition f1model(...) given below
      call hoppetpreEvolve(asQ0,Q0,nloop,xmu,Qmodel)
      call hoppetCachedEvolve(f1model)
      ! alternatively if you need to repeat the same evolution on very
      ! many pdf sets, used a cached evolution (set up once, use many times
      ! and gain a factor 3-4 in speed.
      !call hoppetPreEvolve(asQ0,Q0,nloop,xmu,Q0)
      !call hoppetCachedEvolve(evolvePDF) 
      write(6,'(a)') "Evolution done!"

      ! print out some results
      
      if(ipdf.eq.1) then
        open(unit=1,file='xf1_grid.out') 
      else if(ipdf.eq.2) then
        open(unit=1,file='xf1Tp1_grid.out')
      else if(ipdf.eq.3) then
        open(unit=1,file='xg1_grid.out')
      else if(ipdf.eq.4) then
        open(unit=1,file='xg1T1_grid.out')
      else if(ipdf.eq.5) then
        open(unit=1,file='xh1_grid.out')
      else if(ipdf.eq.6) then
        open(unit=1,file='xh1p1_grid.out')
      else if(ipdf.eq.7) then
        open(unit=1,file='xh1Lp1_grid.out')
      else if(ipdf.eq.8) then
        open(unit=1,file='xh1Tp1_grid.out')
      else if(ipdf.eq.9) then
        open(unit=1,file='xH1p1_grid.out')
      end if
 
*      write(1,'(1p,2a8,11a11)') "#      Q", "x",
*     $     "bbar","cbar","sbar","ubar","dbar","gluon",
*     $     "u","d","s","c","b"

      if(ipdf.le.8) then
      write(6,'(1p,2a11,11a12)') "#      Q", "x",
     $     "bbar","cbar","sbar","ubar","dbar","gluon",
     $     "d","u","s","c","b"
      else if(ipdf.eq.9) then
      write(6,'(1p,2a11,11a12)') "#      Q", "z",
     $     "bbar","cbar","sbar","ubar","dbar","gluon",
     $     "unfav","fav","s","c","b"
      endif


      Qin(1)=Qmin
      Qin(2)=cthre
      Qin(3)=bthre

      Qfin(1)=cthre
      Qfin(2)=bthre
      Qfin(3)=Qmax

C-- new style grid

      do ithre=1,3

         QQmin=Qin(ithre)
         QQmax=Qfin(ithre)

      do iQQ = 1, nQQstep
         !use this for linear grid



         !the following if statement is just to include Qmax in the grid
         if(QQmax.ne.Qmax) then 
         !use this for log grid
         QQi=QQmin*dexp(dble(iQQ-1)*dlog(QQmax/QQmin)/dble(nQQstep))
         else 
         QQi=QQmin*dexp(dble(iQQ-1)*dlog(QQmax/QQmin)/dble(nQQstep-1))
         end if

         xi=xmin-xstep

         do ix = 1, nxstep
         
            xi=xi+xstep

            call hoppetEval(xi,QQi,ourpdf)
            write(6,'(1p,2e11.4,11e12.4)') QQi, xi,
     $        ourpdf(-5), 
     $        ourpdf(-4), 
     $        ourpdf(-3), 
     $        ourpdf(-2), 
     $        ourpdf(-1), 
     $        ourpdf(0),
     $        ourpdf(1), 
     $        ourpdf(2), 
     $        ourpdf(3), 
     $        ourpdf(4),
     $        ourpdf(5)

            write(1,'(1p,2e11.4,11e12.4)') QQi, xi,
     $        ourpdf(-5), 
     $        ourpdf(-4), 
     $        ourpdf(-3), 
     $        ourpdf(-2), 
     $        ourpdf(-1), 
     $        ourpdf(0),
     $        ourpdf(1), 
     $        ourpdf(2), 
     $        ourpdf(3), 
     $        ourpdf(4),
     $        ourpdf(5)

         end do
      end do
      end do

C- old style grid
*
*      QQmin = dsqrt(0.3d0)  
*      QQmax = 10d0
*      nQQstep= 50
*      QQstep = (QQmax-QQmin)/dble(nQQstep-1)
*
*      xmin=0.01
*      xmax=0.9
*      nxstep= 50
*      xstep=(xmax- xmin)/dble(nxstep-1)
*
*
*       !this is just to be sure to include also the first threshold
**       QQi=dsqrt(0.3d0)
*       QQi=1.414213563d0
*
*         xi=xmin-xstep
*
*         do ix = 1, nxstep
*         
*            xi=xi+xstep
*         
*            call hoppetEval(xi,QQi,ourpdf)
*            write(6,'(1p,2e11.4,11e12.4)') QQi, xi,
*     $        ourpdf(-5), 
*     $        ourpdf(-4), 
*     $        ourpdf(-3), 
*     $        ourpdf(-2), 
*     $        ourpdf(-1), 
*     $        ourpdf(0),
*     $        ourpdf(1), 
*     $        ourpdf(2), 
*     $        ourpdf(3), 
*     $        ourpdf(4),
*     $        ourpdf(5)
*
*            write(1,'(1p,2e11.4,11e12.4)') QQi, xi,
*     $        ourpdf(-5), 
*     $        ourpdf(-4), 
*     $        ourpdf(-3), 
*     $        ourpdf(-2), 
*     $        ourpdf(-1), 
*     $        ourpdf(0),
*     $        ourpdf(1), 
*     $        ourpdf(2), 
*     $        ourpdf(3), 
*     $        ourpdf(4),
*     $        ourpdf(5)
*
*         end do
*
*
*      QQi=QQmin-QQstep
*
*      do iQQ = 1, nQQstep
*         !use this for linear grid
*         QQi=QQi+QQstep
*
*         !use this for log grid
*         QQi=QQmin*dexp((iQQ-1)*dlog(QQmax/QQmin)/dble(nQQstep-1))
*
*
*         xi=xmin-xstep
*
*         do ix = 1, nxstep
*         
*            xi=xi+xstep
*
*            call hoppetEval(xi,QQi,ourpdf)
*            write(6,'(1p,2e11.4,11e12.4)') QQi, xi,
*     $        ourpdf(-5), 
*     $        ourpdf(-4), 
*     $        ourpdf(-3), 
*     $        ourpdf(-2), 
*     $        ourpdf(-1), 
*     $        ourpdf(0),
*     $        ourpdf(1), 
*     $        ourpdf(2), 
*     $        ourpdf(3), 
*     $        ourpdf(4),
*     $        ourpdf(5)
*
*            write(1,'(1p,2e11.4,11e12.4)') QQi, xi,
*     $        ourpdf(-5), 
*     $        ourpdf(-4), 
*     $        ourpdf(-3), 
*     $        ourpdf(-2), 
*     $        ourpdf(-1), 
*     $        ourpdf(0),
*     $        ourpdf(1), 
*     $        ourpdf(2), 
*     $        ourpdf(3), 
*     $        ourpdf(4),
*     $        ourpdf(5)
*
*         end do
*      end do
*
*       !this is just to be sure to include also the second threshold
**       QQi=10d0
*       QQi=4.5d0
*
*         xi=xmin-xstep
*
*         do ix = 1, nxstep
*         
*            xi=xi+xstep
*         
*            call hoppetEval(xi,QQi,ourpdf)
*            write(6,'(1p,2e11.4,11e12.4)') QQi, xi,
*     $        ourpdf(-5), 
*     $        ourpdf(-4), 
*     $        ourpdf(-3), 
*     $        ourpdf(-2), 
*     $        ourpdf(-1), 
*     $        ourpdf(0),
*     $        ourpdf(1), 
*     $        ourpdf(2), 
*     $        ourpdf(3), 
*     $        ourpdf(4),
*     $        ourpdf(5)
*
*            write(1,'(1p,2e11.4,11e12.4)') QQi, xi,
*     $        ourpdf(-5), 
*     $        ourpdf(-4), 
*     $        ourpdf(-3), 
*     $        ourpdf(-2), 
*     $        ourpdf(-1), 
*     $        ourpdf(0),
*     $        ourpdf(1), 
*     $        ourpdf(2), 
*     $        ourpdf(3), 
*     $        ourpdf(4),
*     $        ourpdf(5)
*
*         end do
*      print *, hoppetalphas(Qmodel)

      end


      !--------------------------------------------------------------
      ! Model function from Phys.Rev.D78:074010,2008, 
      ! arXiv:0807.0323 [hep-ph].
      !
      ! The subroutine returns x*q(x) in the array pdf(-6:6), containing
      ! flavours (tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t)
      !
      ! The variable Q is set by the calling routine to be equal to the
      ! scale at which the PDF is being requested. Below it is unused
      ! because the initial condition is here provided only for scale
      ! Q0=sqrt(0.3)GeV.
      !
      subroutine f1model(x,Q,pdf)
      double precision x,Q,pdf(-6:6)
      double precision uv, dv
      double precision ubar, dbar
      !---------------------
      
      double precision mq,Lams,Ms,Lamv,Mv,Lamvv,Mvv
      double precision gs, gv, gvv

      external f1sint
      external f1vint
      external f1vvint

*      external H1p1fav
*      external H1p1favint
      double precision z
      common/collins/z

      common/pdftype/ipdf
      common/f1int/mq,Lams,Ms,Lamv,Mv,Lamvv,Mvv

      data mq,Lams,Ms/0.3, 0.609, 0.822/
      data Lamv,Mv,Lamvv,Mvv/0.716, 1.492, 0.376, 0.890/
      data gs, gv, gvv/0.847, 1.061, 0.880/

      if(ipdf.eq.1) then
      uv = x*(gs**2*f1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*f1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*f1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.2) then
      uv = x*(gs**2*f1Tp1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*f1Tp1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*f1Tp1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.3) then
      uv = x*(gs**2*g1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*g1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*g1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.4) then
      uv = x*(gs**2*g1T1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*g1T1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*g1T1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.5) then
      uv = x*(gs**2*h1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*h1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*h1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.6) then
      uv = x*(gs**2*h1p1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*h1p1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*h1p1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.7) then
      uv = x*(gs**2*h1Lp1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*h1Lp1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*h1Lp1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.8) then
      uv = x*(gs**2*h1Tp1sunnorm(x,mq,Lams,Ms)/
     +             dgauss(f1sint,0d0,1d0,1.d-2)
     -            +gv**2*h1Tp1vunnorm(x,mq,Lamv,Mv)/
     +          dgauss(f1vint,0d0,1d0,1.d-2))

      dv = x*gvv**2*h1Tp1vunnorm(x,mq,Lamvv,Mvv)/
     +            dgauss(f1vvint,0d0,1d0,1.d-2)
      else if(ipdf.eq.9) then
      z=x
      uv = x*dgauss(H1p1fav,0d0,1d0,1.d-2)
     &     +x*dgauss(H1p1favint,0d0,1d0,1.d-2)
      dv = -uv
      end if

      dbar = 0
      ubar = 0

      pdf(0) = 0
      pdf(-3) = 0
      pdf( 3) = 0
      pdf( 2) = uv
      pdf(-2) = 0
      pdf( 1) = dv
      pdf(-1) = 0

      pdf( 4) = 0
      pdf( 5) = 0
      pdf( 6) = 0
      pdf(-4) = 0
      pdf(-5) = 0
      pdf(-6) = 0
      end 

C----------------------------------------------------------------------
      real*8 function f1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      f1sunnorm=((1 - x)**3*(2*(mq + M*x)**2 + Ls2))/
     -  (96.*Pi**2*Ls2**3)

      return
      end

C----------------------------------------------------------------------
      real*8 function f1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      f1vunnorm=((1 - x)*(2*(mq + M*x)**2*(1 - x)**2 + 
     -      (1 + x**2)*Lv2))/
     -  (96.*Pi**2*Lv2**3)

      return
      end

C----------------------------------------------------------------------
      real*8 function f1Tp1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      common/alphas/Qmodel

*      asQmodel = 0.3

      asQmodel = hoppetalphas(Qmodel)

      print *, asQmodel

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      f1Tp1sunnorm=-4*Pi*4/3*asQmodel*((mq + M*x)*(1-x)**3)/
     -  (256.*Pi**3*Ls2**2*M)

      return
      end

C---------------------------------------------------------------------- 
      real*8 function f1Tp1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      common/alphas/Qmodel

*      asQmodel = 0.3

      asQmodel = hoppetalphas(Qmodel)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      f1Tp1vunnorm=4*Pi*4/3*asQmodel*(x*(mq + M*x)*(1-x)**2)/
     -  (256.*Pi**3*Lv2**2*M)

      return
      end

C----------------------------------------------------------------------
      real*8 function g1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      g1sunnorm=((1 - x)**3*(2*(mq + M*x)**2 - Ls2))/
     -  (96.*Pi**2*Ls2**3)

      return
      end

C---------------------------------------------------------------------- 
      real*8 function g1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      g1vunnorm=-((1 - x)*(2*(mq + M*x)**2*(1 - x)**2 - 
     -      (1 + x**2)*Lv2))/
     -  (96.*Pi**2*Lv2**3)

      return
      end


C----------------------------------------------------------------------
      real*8 function g1T1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      g1T1sunnorm= M*(mq + M*x)*(1 - x)**3/
     -  (2.*M**2*48.*Pi**2*Ls2**2)

      return
      end

C----------------------------------------------------------------------
      real*8 function g1T1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      g1T1vunnorm= M*(mq + M*x)*x*(1 - x)**2/
     -  (2.*M**2*48.*Pi**2*Lv2**2)

      return
      end

C----------------------------------------------------------------------
      real*8 function h1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      h1sunnorm=((1 - x)**3*(mq + M*x)**2)/
     -  (48.*Pi**2*Ls2**3)

      return
      end

C----------------------------------------------------------------------
      real*8 function h1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      h1vunnorm=- x*(1 - x)/
     -  (48.*Pi**2*Lv2**2)

      return
      end


C----------------------------------------------------------------------
      real*8 function h1p1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      common/alphas/Qmodel

      data Pi/3.141592654d0/
      data M/0.938d0/

*      asQmodel=0.3

      asQmodel = hoppetalphas(Qmodel)

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      h1p1sunnorm=-4*Pi*4/3*asQmodel*((mq + M*x)*(1-x)**3)/
     -  (256.*Pi**3*Ls2**2*M)

      return
      end

C----------------------------------------------------------------------
      real*8 function h1p1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      common/alphas/Qmodel

      data Pi/3.141592654d0/
      data M/0.938d0/

*      asQmodel=0.3

      asQmodel = hoppetalphas(Qmodel)

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      h1p1vunnorm=-4*Pi*4/3*asQmodel*((mq + M*x)*(1-x)**2)/
     -  (256.*Pi**3*Lv2**2*M)

      return
      end

C----------------------------------------------------------------------
      real*8 function h1Lp1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      h1Lp1sunnorm=- M*(mq + M*x)*(1 - x)**3/
     -  (2.*M**2*48.*Pi**2*Ls2**2)

      return
      end

C---------------------------------------------------------------------- 
      real*8 function h1Lp1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      h1Lp1vunnorm= M*(mq + M*x)*(1 - x)**2/
     -  (2.*M**2*48.*Pi**2*Lv2**2)

      return
      end

C----------------------------------------------------------------------
      real*8 function h1Tp1sunnorm(x,mq,Lams,Ms)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Ls2=Lams**2*(1 - x) + Ms**2*x - M**2*(1 - x)*x

      h1Tp1sunnorm=- M**2*(1 - x)**3/
     -  (2.*M**2*48.*Pi**2*Ls2**2)

      return
      end

C---------------------------------------------------------------------- 
      real*8 function h1Tp1vunnorm(x,mq,Lamv,Mv)
      implicit real*8 (a-h,k-z)

      data Pi/3.141592654d0/
      data M/0.938d0/

      Lv2=Lamv**2*(1 - x) + Mv**2*x - M**2*(1 - x)*x

      h1Tp1vunnorm= 0d0

      return
      end


C----------------------------------------------------------------------
C   The following functions are needed for some integrations
C----------------------------------------------------------------------

      real*8 function f1sint(x)
      implicit real*8 (a-h,k-z)
      
      common/f1int/mq,Lams,Ms,Lamv,Mv,Lamvv,Mvv

      f1sint=f1sunnorm(x,mq,Lams,Ms)

      return
      end


      real*8 function f1vint(x)
      implicit real*8 (a-h,k-z)
      
      common/f1int/mq,Lams,Ms,Lamv,Mv,Lamvv,Mvv

      f1vint=f1vunnorm(x,mq,Lamv,Mv)

      return
      end


      real*8 function f1vvint(x)
      implicit real*8 (a-h,k-z)
      
      common/f1int/mq,Lams,Ms,Lamv,Mv,Lamvv,Mvv

      f1vvint=f1vunnorm(x,mq,Lamvv,Mvv)

      return
      end
