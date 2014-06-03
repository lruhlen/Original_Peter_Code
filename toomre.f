

      program toomre

c..   Computes galaxy encounters using the restricted three-body approximation
c..   Requires Numerical Recipes Routines: RAN2 
c.                                         BSSTEP  (set nmax=2000)
c.                                         PZEXTR  (set nmax=2000)
c.                                         MMID    (set nmax=2000)
c...
c............................................................................

      implicit real*8(a-h,o-z), integer*4(i-n)
      parameter (nb=176,ngal=2,numrings=6)
      parameter(n=6*nb,pi=3.14159265358979323846)
      dimension nring(6)
      dimension rring(6)
      dimension y(n),dydx(n),yscal(n)
      common /path2/ rm(ngal),vcx,vcy,vcz
      data nring /14,20,26,32,38,44/
      data rring /.14,.20,.26,.32,.38,.44/
      external derivs

c..   Timestep accuracy for Bulirsch-Stoer
      Acc=1.e-8
c..   number of phase space dimensions
      nvar=n
c..   number of scattering trials
      nloop=1
c..   initial timestep guess
      htry=0.01
c..   integration timescale 
      x2=7.
c..   maximum timestep number
      nstepmax=5000
c..   printout interval
      printdt=1.
c..   starting time for the integration
      x=0.
c..   initial snapshot number
      nsnap=0

      open(unit=1,file='gal1.pos',status='unknown')
      open(unit=2,file='gal2.pos',status='unknown')

c..   Galaxy encounter initial conditions:

      rm(1)=5.
      rm(2)=1.
      a=3.
      e=0.6
      rinc=pi/6.
      p=pi/4.
      rn=pi/2.
      rl=1.7*pi

c..   euler angles for rotation of angular momentum vector of initial galaxy
      ephi=pi/3.
      ethe=pi/3.
      epsi=0.

c..   transform orbital elements of galaxy pair into cartesian coordinates:
      y(1)=0.
      y(2)=0.
      y(3)=0.
      y(4)=0.
      y(5)=0.
      y(6)=0.

      rmu=rm(1)+rm(2)
      q=a*(1-e)

c..   this routine does the elements to cartesian conversion
      call mco_el2x (rmu,q,e,rinc,p,rn,rl,rx,ry,rz,ru,rv,rw)
      ib=6*(ngal-1)

      y(ib+1)=rx
      y(ib+3)=ry
      y(ib+5)=rz

      y(ib+2)=ru
      y(ib+4)=rv
      y(ib+6)=rw


c..   form the rotation matrix for the given euler angles:
      s1 = sin(ephi)
      c1 = cos(ephi)
      s2 = sin(ethe)
      c2 = cos(ethe)
      s3 = sin(epsi)
      c3 = cos(epsi)
      a11 = c3*c1 - c2*s1*s3
      a12 = -1.0*s3*c1 - c2*s1*c3
      a13 = s2*s1
      a21 = c3*s1 + c2*c1*s3
      a22 = -1.0*s3*s1 + c2*c1*c3
      a23 = -1.0*s2*c1
      a31 = s3*s2
      a32 = c3*s2
      a33 = c2

c..   set up the test particles
      ipart=ngal
      do j=1,numrings
        do i=1,nring(j)

          ipart=ipart+1
          ib=(ipart-1)*6
          thet=(real(i)/real(nring(j)))*(2.*pi)

          xp0=rring(j)*cos(thet)
          yp0=rring(j)*sin(thet)
          zp0=0.0

          y(ib+1)=a11*xp0 + a12*yp0 + a13*zp0
          y(ib+3)=a21*xp0 + a22*yp0 + a23*zp0
          y(ib+5)=a31*xp0 + a32*yp0 + a33*zp0

          vcirc= sqrt(rm(2)/rring(j))

          vx0=-vcirc*sin(thet)
          vy0= vcirc*cos(thet)
          vz0=0.0

          y(ib+2)=a11*vx0 + a12*vy0 + a13*vz0
          y(ib+4)=a21*vx0 + a22*vy0 + a23*vz0
          y(ib+6)=a31*vx0 + a32*vy0 + a33*vz0

        end do
      end do

c..   put them in orbit around the primary
      do i=ngal+1,ipart
        ib=(i-1)*6
        kb=(ngal-1)*6
        do jb=1,6
          y(ib+jb)=y(ib+jb)+y(kb+jb)
        end do
      end do

c..   Compute center of mass velocity for the system and subtract it.
     
      vcx=0.
      vcy=0.
      vcz=0.

      rmtot=0.
      do i=1,ngal
        ib=(i-1)*6
        rmtot=rmtot+rm(i)
        vcx=vcx+rm(i)*y(ib+2)
        vcy=vcy+rm(i)*y(ib+4)
        vcz=vcz+rm(i)*y(ib+6)
      end do

      vcx=vcx/rmtot
      vcy=vcy/rmtot
      vcz=vcz/rmtot

      do i=1,nb
        ib=(i-1)*6
        y(ib+2)=y(ib+2)-vcx
        y(ib+4)=y(ib+4)-vcy
        y(ib+6)=y(ib+6)-vcz
      end do

c..   Compute center of mass position of system and subtract it.
      pcx=0.
      pcy=0.
      pcz=0.

      rmtot=0.
      do i=1,ngal
        ib=(i-1)*6
        rmtot=rmtot+rm(i)
        pcx=pcx+rm(i)*y(ib+1)
        pcy=pcy+rm(i)*y(ib+3)
        pcz=pcz+rm(i)*y(ib+5)
      end do

      pcx=pcx/rmtot
      pcy=pcy/rmtot
      pcz=pcz/rmtot

      do i=1,nb
        ib=(i-1)*6
        y(ib+1)=y(ib+1)-pcx
        y(ib+3)=y(ib+3)-pcy
        y(ib+5)=y(ib+5)-pcz
      end do


c..   Overall integration loop
      iflag=0
2105  continue

        iflag=iflag+1
        call derivs(x,y,dydx)
        do iscale=1,nvar
          yscal(iscale)=abs(y(iscale))+abs(htry*dydx(iscale))+1.e-30
        end do
        call bsstep(y,dydx,nvar,x,htry,acc,yscal,hdid,hnext,derivs)  
        htry=hnext
     
c..     write out high-resolution orbit for the two test particles:
        write(1,*) y(1),y(3),y(5)
        write(2,*) y(7),y(9),y(11)

        if(x.gt.(nsnap*printdt)) then
          nsnap=nsnap+1
          call printout(x,y,nsnap)
        end if

        if(iflag.gt.nstepmax.or.x.gt.x2) then
          goto 2106
        end if

c..   go do another step.
      goto 2105

2106  continue
c..   End of overall integration loop for a particular experiment

      end


      subroutine derivs(x,y,dydx)

      implicit real*8(a-h,o-z), integer*4(i-n)

      parameter (nb=176,ngal=2)
      parameter (eps2=0.0005)
      parameter(n=6*nb,pi=3.14159265358979323846)


      common /path2/ rm(ngal),vcx,vcy,vcz


      real*8 x,y(n),dydx(n)
      real*8 denom(ngal,ngal)
      real*8 denom2(nb,ngal)


      do i=2,n,2
        dydx(i-1)=y(i)
      end do


      do i=1,ngal 
        do j=i+1,ngal 

          if (i.ne.j) then
            jb=(j-1)*6
            ib=(i-1)*6
            denom(i,j)=(y(jb+1)-y(ib+1))*(y(jb+1)-y(ib+1)) +
     +                 (y(jb+3)-y(ib+3))*(y(jb+3)-y(ib+3)) +
     +                 (y(jb+5)-y(ib+5))*(y(jb+5)-y(ib+5)) +
     +                 eps2

            denom(i,j)=denom(i,j)*sqrt(denom(i,j))
            denom(j,i)=denom(i,j)
          end if

        end do
      end do

      do i=ngal +1,nb
        do j=1,ngal 

          if (i.ne.j) then
            jb=(j-1)*6
            ib=(i-1)*6
            denom2(i,j)=(y(jb+1)-y(ib+1))*(y(jb+1)-y(ib+1)) +
     +                  (y(jb+3)-y(ib+3))*(y(jb+3)-y(ib+3)) +
     +                  (y(jb+5)-y(ib+5))*(y(jb+5)-y(ib+5)) +
     +                  eps2

            denom2(i,j)=denom2(i,j)*sqrt(denom2(i,j))
          end if

        end do
      end do

2104  continue

c..   self-gravitating particles
      do i=1,ngal 
        ib=(i-1)*6

        do ic=1,3

         dydx(ib+(2*ic))=0.

c..       contributions from the other massive particles:
          do j=1,ngal 
            jb=(j-1)*6
            if( i.ne.j) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic)) -
     +        rm(j)*(y(ib+2*ic-1)-y(jb+2*ic-1))/denom(i,j)
            end if
          end do

        end do
      end do

c..   test particles
      do i=ngal +1,nb
        ib=(i-1)*6

        do ic=1,3

          dydx(ib+(2*ic))=0.

c..       contributions from the massive particles:
          do j=1,ngal 
            jb=(j-1)*6
            dydx(ib+(2*ic))=dydx(ib+(2*ic)) -
     +      rm(j)*(y(ib+2*ic-1)-y(jb+2*ic-1))/denom2(i,j)
          end do

        end do
      end do






 
      return
      end



c..................................................................
c
      subroutine printout(x,y,nsnap)
c
c..................................................................

       implicit real*8(a-h,o-z), integer*4(i-n)

       parameter (nb=176,ngal=2)
       parameter(n=6*nb,pi=3.14159265358979323846)


      common /path2/ rm(ngal),vcx,vcy,vcz
      real*8 y(6*nb)

      character*3 sstring
      character*4 filestem
      character*10 nstring
      character*7 target_out

      data nstring/'0123456789'/
      save
      sstring(1:1)=nstring(1+nsnap/100:1+nsnap/100)
      istring=1+mod(nsnap,100)/10
      sstring(2:2)=nstring(istring:istring)
      istring=1+mod(nsnap,10)
      sstring(3:3)=nstring(istring:istring)
      filestem='snap'
      target_out=filestem(1:4)//sstring(1:3)
      open(unit=10,file=target_out,status='unknown')

      do i=1,nb
        jb=(i-1)*6
        write(10,2104) y(jb+1),y(jb+3),y(jb+5)
      end do
2104  format(4e13.5)

      close(10)

      return
      end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_X2EL.FOR    (ErikSoft  6 May 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates Keplerian orbital elements given relative coordinates and
c velocities, and MU = G times the sum of the masses.
c
c The elements are: q = perihelion distance
c                   e = eccentricity
c                   i = inclination
c                   p = longitude of perihelion (NOT argument of perihelion!!)
c                   n = longitude of ascending node
c                   l = mean anomaly (or mean longitude if e < 1.e-8)
c
c------------------------------------------------------------------------------
c
      subroutine mco_x2el (mu,x,y,z,u,v,w,q,e,i,p,n,l)
c
      implicit none
      integer NMAX, CMAX, NMESS
      real*8 HUGE

      parameter (NMAX = 2000)
      parameter (CMAX = 50)
      parameter (NMESS = 200)
      parameter (HUGE = 9.9d29)

c Constants:
c
c DR = conversion factor from degrees to radians
c K2 = Gaussian gravitational constant squared
c AU = astronomical unit in cm
c MSUN = mass of the Sun in g
c
      real*8 PI,TWOPI,PIBY2,DR,K2,AU,MSUN
c
      parameter (PI = 3.141592653589793d0)
      parameter (TWOPI = PI * 2.d0)
      parameter (PIBY2 = PI * .5d0)
      parameter (DR = PI / 180.d0)
      parameter (K2 = 2.959122082855911d-4)
      parameter (AU = 1.4959787e13)
      parameter (MSUN = 1.9891e33)


c
c Input/Output
      real*8 mu,q,e,i,p,n,l,x,y,z,u,v,w
c
c Local
      real*8 hx,hy,hz,h2,h,v2,r,rv,s,true
      real*8 ci,to,temp,tmp2,bige,f,cf,ce
c
c------------------------------------------------------------------------------
c
      hx = y * w  -  z * v
      hy = z * u  -  x * w
      hz = x * v  -  y * u
      h2 = hx*hx + hy*hy + hz*hz
      v2 = u * u  +  v * v  +  w * w
      rv = x * u  +  y * v  +  z * w
      r = sqrt(x*x + y*y + z*z)
      h = sqrt(h2)
      s = h2 / mu
c
c Inclination and node
      ci = hz / h
      if (abs(ci).lt.1) then
        i = acos (ci)
        n = atan2 (hx,-hy)
        if (n.lt.0) n = n + TWOPI
      else
        if (ci.gt.0) i = 0.d0
        if (ci.lt.0) i = PI
        n = 0.d0
      end if
c
c Eccentricity and perihelion distance
      temp = 1.d0 + s*(v2/mu - 2.d0/r)
      if (temp.le.0) then
        e = 0.d0
      else
        e = sqrt (temp)
      end if
      q = s / (1.d0 + e)
c
c True longitude
      if (hy.ne.0) then
        to = -hx/hy
        temp = (1.d0 - ci) * to
        tmp2 = to * to
        true = atan2((y*(1.d0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
      else
        true = atan2(y * ci, x)
      end if
      if (ci.lt.0) true = true + PI
c
      if (e.lt.1.d-8) then
        p = 0.d0
        l = true
      else
        ce = (v2*r - mu) / (e*mu)
c
c Mean anomaly for ellipse
        if (e.lt.1) then
          if (abs(ce).gt.1) ce = sign(1.d0,ce)
          bige = acos(ce)
          if (rv.lt.0) bige = TWOPI - bige
          l = bige - e*sin(bige)
        else
c
c Mean anomaly for hyperbola
          if (ce.lt.1) ce = 1.d0
          bige = log( ce + sqrt(ce*ce-1.d0) )
          if (rv.lt.0) bige = TWOPI - bige
          l = e*sinh(bige) - bige
        end if
c
c Longitude of perihelion
        cf = (s - r) / (e*r)
        if (abs(cf).gt.1) cf = sign(1.d0,cf)
        f = acos(cf)
        if (rv.lt.0) f = TWOPI - f
        p = true - f
        p = mod (p + TWOPI + TWOPI, TWOPI)
      end if
c
      if (l.lt.0) l = l + TWOPI
      if (l.gt.TWOPI) l = mod (l, TWOPI)
c
c------------------------------------------------------------------------------
c
      return
      end



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_KEP.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Solves Kepler's equation for eccentricities less than one.
c Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
c
c  e = eccentricity
c  l = mean anomaly      (radians)
c  u = eccentric anomaly (   "   )
c
c------------------------------------------------------------------------------
c
      function mco_kep (e,oldl)
      implicit none
c
c Input/Outout
      real*8 oldl,e,mco_kep
c
c Local
      real*8 l,pi,twopi,piby2,u1,u2,ome,sign
      real*8 x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
      real*8 p,q,p2,ss,cc
      logical flag,big,bigg
c
c------------------------------------------------------------------------------
c
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
      piby2 = .5d0 * pi
c
c Reduce mean anomaly to lie in the range 0 < l < pi
      if (oldl.ge.0) then
        l = mod(oldl, twopi)
      else
        l = mod(oldl, twopi) + twopi
      end if
      sign = 1.d0
      if (l.gt.pi) then
        l = twopi - l
        sign = -1.d0
      end if
c
      ome = 1.d0 - e
c
      if (l.ge..45d0.or.e.lt..55d0) then
c
c Regions A,B or C in Nijenhuis
c -----------------------------
c
c Rough starting value for eccentric anomaly
        if (l.lt.ome) then
          u1 = ome
        else
          if (l.gt.(pi-1.d0-e)) then
            u1 = (l+e*pi)/(1.d0+e)
          else
            u1 = l + e
          end if
        end if
c
c Improved value using Halley's method
        flag = u1.gt.piby2
        if (flag) then
          x = pi - u1
        else
          x = u1
        end if
        x2 = x*x
        sn = x*(1.d0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.d0 + x2*(-.49815 + x2*.03805)
        if (flag) dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.d0 - e*dsn
        u2 = u1 - f0/(f1 - .5d0*f0*f2/f1)
      else
c
c Region D in Nijenhuis
c ---------------------
c
c Rough starting value for eccentric anomaly
        z1 = 4.d0*e + .5d0
        p = ome / z1
        q = .5d0 * l / z1
        p2 = p*p
        z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.d0*q / ( z2 + p + p2/z2 )
c
c Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3)
        u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 )
      end if
c
c Accurate value using 3rd-order version of Newton's method
c N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
c
c First get accurate values for u2 - sin(u2) and 1 - cos(u2)
      bigg = (u2.gt.piby2)
      if (bigg) then
        z3 = pi - u2
      else
        z3 = u2
      end if
c
      big = (z3.gt.(.5d0*piby2))
      if (big) then
        x = piby2 - z3
      else
        x = z3
      end if
c
      x2 = x*x
      ss = 1.d0
      cc = 1.d0
c
      ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. -
     %   x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
      cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. -
     %   x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. -
     %   x2/306.))))))))
c
      if (big) then
        z1 = cc + z3 - 1.d0
        z2 = ss + z3 + 1.d0 - piby2
      else
        z1 = ss
        z2 = cc
      end if
c
      if (bigg) then
        z1 = 2.d0*u2 + z1 - pi
        z2 = 2.d0 - z2
      end if
c
      f0 = l - u2*ome - e*z1
      f1 = ome + e*z2
      f2 = .5d0*e*(u2-z1)
      f3 = e/6.d0*(1.d0-z2)
      z1 = f0/f1
      z2 = f0/(f2*z1+f1)
      mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
c
c------------------------------------------------------------------------------
c
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_SINE.FOR    (ErikSoft  17 April 1997)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates sin and cos of an angle X (in radians).
c
c------------------------------------------------------------------------------
c
      subroutine mco_sine (x,sx,cx)
c
      implicit none
c
c Input/Output
      real*8 x,sx,cx
c
c Local
      real*8 pi,twopi
c
c------------------------------------------------------------------------------
c
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
c
      if (x.gt.0) then
        x = mod(x,twopi)
      else
        x = mod(x,twopi) + twopi
      end if
c
      cx = cos(x)
c
      if (x.gt.pi) then
        sx = -sqrt(1.d0 - cx*cx)
      else
        sx =  sqrt(1.d0 - cx*cx)
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_SINH.FOR    (ErikSoft  12 June 1998)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Calculates sinh and cosh of an angle X (in radians)
c
c------------------------------------------------------------------------------
c
      subroutine mco_sinh (x,sx,cx)
c
      implicit none
c
c Input/Output
      real*8 x,sx,cx
c
c------------------------------------------------------------------------------
c
      sx = sinh(x)
      cx = sqrt (1.d0 + sx*sx)
c
c------------------------------------------------------------------------------
c
      return
      end
***********************************************************************
c                    ORBEL_FGET.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_fget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
*           Cel. Mech. ".  Quartic convergence from Danby's book.
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_fget(e,capn)

        implicit NONE    

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
       parameter  (NTPMAX = 1001) ! max number of test particles

c...   Size of the test particle integer status flag
        integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)
        integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

c...   Size of the test particle integer status flag
        integer NSTATR
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
        real*8 DANBYAC , DANBYB
        parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
        integer NLAG1, NLAG2
        parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
        real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
        real*8 PI,TWOPI,PIBY2,DEGRAD
        parameter (PI = 3.14159265358979D0)
        parameter (TWOPI = 2.0D0 * PI)
        parameter (PIBY2 = PI/2.0D0)
        parameter (DEGRAD = 180.0D0 / PI)


c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer i,IMAX
	real*8 tmp,x,shx,chx
	real*8 esh,ech,f,fp,fpp,fppp,dx
	PARAMETER (IMAX = 10)

c----
c...  Executable code 

c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. 

c  begin with a guess proposed by Danby	
	if( capn .lt. 0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0
	   x = -log(tmp)
	else
	   tmp = +2.d0*capn/e + 1.8d0
	   x = log( tmp)
	endif

	orbel_fget = x

	do i = 1,IMAX
	  call orbel_schget(x,shx,chx)
	  esh = e*shx
	  ech = e*chx
	  f = esh - x - capn
c	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_fget = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) RETURN
	  x = orbel_fget
	enddo	

	write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	return
	end   ! orbel_fget
c------------------------------------------------------------------
***********************************************************************
c                    ORBEL_FLON.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_flon ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Uses power series for N in terms of F and Newton,s method
*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26, 1992.
*     REVISIONS: 
***********************************************************************

	real*8 function orbel_flon(e,capn)

        implicit NONE    

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
       parameter  (NTPMAX = 1001) ! max number of test particles

c...   Size of the test particle integer status flag
        integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)
        integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

c...   Size of the test particle integer status flag
        integer NSTATR
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
        real*8 DANBYAC , DANBYB
        parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
        integer NLAG1, NLAG2
        parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
        real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
        real*8 PI,TWOPI,PIBY2,DEGRAD
        parameter (PI = 3.14159265358979D0)
        parameter (TWOPI = 2.0D0 * PI)
        parameter (PIBY2 = PI/2.0D0)
        parameter (DEGRAD = 180.0D0 / PI)


c...  Inputs Only: 
	real*8 e,capn

c...  Internals:
	integer iflag,i,IMAX
	real*8 a,b,sq,biga,bigb
	real*8 x,x2
	real*8 f,fp,dx
	real*8 diff
	real*8 a0,a1,a3,a5,a7,a9,a11
	real*8 b1,b3,b5,b7,b9,b11
	PARAMETER (IMAX = 10)
	PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
	PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
	PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
	PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

c----
c...  Executable code 


c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. Only good for smallish CAPN 

	iflag = 0
	if( capn .lt. 0.d0) then
	   iflag = 1
	   capn = -capn
	endif

	a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
	a0 = -6227020800.d0*capn/e
	b1 = a1

c  Set iflag nonzero if capn < 0., in which case solve for -capn
c  and change the sign of the final answer for F.
c  Begin with a reasonable guess based on solving the cubic for small F	


	a = 6.d0*(e-1.d0)/e
	b = -6.d0*capn/e
	sq = sqrt(0.25*b*b +a*a*a/27.d0)
	biga = (-0.5*b + sq)**0.3333333333333333d0
	bigb = -(+0.5*b + sq)**0.3333333333333333d0
	x = biga + bigb
c	write(6,*) 'cubic = ',x**3 +a*x +b
	orbel_flon = x
c If capn is tiny (or zero) no need to go further than cubic even for
c e =1.
	if( capn .lt. TINY) go to 100

	do i = 1,IMAX
	  x2 = x*x
	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
	  dx = -f/fp
c	  write(6,*) 'i,dx,x,f : '
c	  write(6,432) i,dx,x,f
432	  format(1x,i3,3(2x,1p1e22.15))
	  orbel_flon = x + dx
c   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) go to 100
	  x = orbel_flon
	enddo	

c Abnormal return here - we've gone thru the loop 
c IMAX times without convergence
	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif
	write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	  diff = e*sinh(orbel_flon) - orbel_flon - capn
	  write(6,*) 'N, F, ecc*sinh(F) - F - N : '
	  write(6,*) capn,orbel_flon,diff
	return

c  Normal return here, but check if capn was originally negative
100	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif

	return
	end     ! orbel_flon
c------------------------------------------------------------------
***********************************************************************
c	                  ORBEL_SCGET.F
***********************************************************************
*     PURPOSE:  Given an angle, efficiently compute sin and cos.
*
*        Input:
*             angle ==> angle in radians (real scalar)
*        
*        Output:
*             sx    ==>  sin(angle)  (real scalar)
*             cx    ==>  cos(angle)  (real scalar)
*
*     ALGORITHM: Obvious from the code 
*     REMARKS: The HP 700 series won't return correct answers for sin
*       and cos if the angle is bigger than 3e7. We first reduce it
*       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
*       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 6, 1992.
*     REVISIONS: 
***********************************************************************

	subroutine orbel_scget(angle,sx,cx)

        implicit NONE    

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
       parameter  (NTPMAX = 1001) ! max number of test particles

c...   Size of the test particle integer status flag
        integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)
        integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

c...   Size of the test particle integer status flag
        integer NSTATR
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
        real*8 DANBYAC , DANBYB
        parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
        integer NLAG1, NLAG2
        parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
        real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
        real*8 PI,TWOPI,PIBY2,DEGRAD
        parameter (PI = 3.14159265358979D0)
        parameter (TWOPI = 2.0D0 * PI)
        parameter (PIBY2 = PI/2.0D0)
        parameter (DEGRAD = 180.0D0 / PI)


c...  Inputs Only: 
        real*8 angle

c...  Output:
	real*8 sx,cx

c... Internals:
	integer nper
	real*8 x
	real*8 PI3BY2
	parameter(PI3BY2 = 1.5d0*PI)

c----
c...  Executable code 

        nper = angle/TWOPI
	x = angle - nper*TWOPI
	if(x.lt.0.d0) then
           x = x + TWOPI
        endif
	sx = sin(x)
	cx= sqrt(1.d0 - sx*sx)
	if( (x .gt. PIBY2) .and. (x .lt.PI3BY2)) then
           cx = -cx
        endif

	return
	end   ! orbel_scget
c-------------------------------------------------------------------
***********************************************************************
c	                  ORBEL_SCHGET.F
***********************************************************************
*     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
*
*        Input:
*             angle ==> angle in radians (real scalar)
*        
*        Output:
*             shx    ==>  sinh(angle)  (real scalar)
*             chx    ==>  cosh(angle)  (real scalar)
*
*     ALGORITHM: Obvious from the code 
*     REMARKS: Based on the routine SCGET for sine's and cosine's.
*       We use the sqrt rather than cosh (it's faster)
*       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
*       OR OVERFLOWS WILL OCCUR!
*     AUTHOR:  M. Duncan.
*     DATE WRITTEN:  May 6, 1992.
*     REVISIONS: 
***********************************************************************

	subroutine orbel_schget(angle,shx,chx)

        implicit NONE    

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
       parameter  (NTPMAX = 1001) ! max number of test particles

c...   Size of the test particle integer status flag
        integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)
        integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

c...   Size of the test particle integer status flag
        integer NSTATR
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
        real*8 DANBYAC , DANBYB
        parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
        integer NLAG1, NLAG2
        parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
        real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
        real*8 PI,TWOPI,PIBY2,DEGRAD
        parameter (PI = 3.14159265358979D0)
        parameter (TWOPI = 2.0D0 * PI)
        parameter (PIBY2 = PI/2.0D0)
        parameter (DEGRAD = 180.0D0 / PI)


c...  Inputs Only: 
        real*8 angle

c...  Output:
	real*8 shx,chx

c----
c...  Executable code 

	shx = sinh(angle)
	chx= sqrt(1.d0 + shx*shx)

	return
	end   ! orbel_schget
c---------------------------------------------------------------------
***********************************************************************
c                    ORBEL_ZGET.F
***********************************************************************
*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
*          given Q (Fitz. notation.)
*
*             Input:
*                           q ==>  parabola mean anomaly. (real scalar)
*             Returns:
*                  orbel_zget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
*     REMARKS: For a parabola we can solve analytically.
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: May 27 - corrected it for negative Q and use power
*	      series for small Q.
***********************************************************************

	real*8 function orbel_zget(q)

            implicit NONE    

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
       parameter  (NTPMAX = 1001) ! max number of test particles

c...   Size of the test particle integer status flag
        integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)
        integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

c...   Size of the test particle integer status flag
        integer NSTATR
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
        real*8 DANBYAC , DANBYB
        parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
        integer NLAG1, NLAG2
        parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
        real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
        real*8 PI,TWOPI,PIBY2,DEGRAD
        parameter (PI = 3.14159265358979D0)
        parameter (TWOPI = 2.0D0 * PI)
        parameter (PIBY2 = PI/2.0D0)
        parameter (DEGRAD = 180.0D0 / PI)


c...  Inputs Only: 
	real*8 q

c...  Internals:
	integer iflag
	real*8 x,tmp

c----
c...  Executable code 

	iflag = 0
	if(q.lt.0.d0) then
	  iflag = 1
	  q = -q
	endif

	if (q.lt.1.d-3) then
	   orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
	else
	   x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
	   tmp = x**(1.d0/3.d0)
	   orbel_zget = tmp - 1.d0/tmp
	endif

	if(iflag .eq.1) then
           orbel_zget = -orbel_zget
	   q = -q
	endif
	
	return
	end    ! orbel_zget
c----------------------------------------------------------------------



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_EL2X.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates Cartesian coordinates and velocities given Keplerian orbital
c elements (for elliptical, parabolic or hyperbolic orbits).
c
c Based on a routine from Levison and Duncan's SWIFT integrator.
c
c  mu = grav const * (central + secondary mass)
c  q = perihelion distance
c  e = eccentricity
c  i = inclination                 )
c  p = longitude of perihelion !!! )   in
c  n = longitude of ascending node ) radians
c  l = mean anomaly                )
c
c  x,y,z = Cartesian positions  ( units the same as a )
c  u,v,w =     "     velocities ( units the same as sqrt(mu/a) )
c
c------------------------------------------------------------------------------
c
      subroutine mco_el2x (mu,q,e,i,p,n,l,x,y,z,u,v,w)
c
      implicit none
      integer NMAX, CMAX, NMESS
      real*8 HUGE

      parameter (NMAX = 2000)
      parameter (CMAX = 50)
      parameter (NMESS = 200)
      parameter (HUGE = 9.9d29)

c Constants:
c
c DR = conversion factor from degrees to radians
c K2 = Gaussian gravitational constant squared
c AU = astronomical unit in cm
c MSUN = mass of the Sun in g
c
      real*8 PI,TWOPI,PIBY2,DR,K2,AU,MSUN
c
      parameter (PI = 3.141592653589793d0)
      parameter (TWOPI = PI * 2.d0)
      parameter (PIBY2 = PI * .5d0)
      parameter (DR = PI / 180.d0)
      parameter (K2 = 2.959122082855911d-4)
      parameter (AU = 1.4959787e13)
      parameter (MSUN = 1.9891e33)


c
c Input/Output
      real*8 mu,q,e,i,p,n,l,x,y,z,u,v,w
c
c Local
      real*8 g,a,ci,si,cn,sn,cg,sg,ce,se,romes,temp
      real*8 z1,z2,z3,z4,d11,d12,d13,d21,d22,d23
      real*8 mco_kep, orbel_fhybrid, orbel_zget
c
c------------------------------------------------------------------------------
c
c Change from longitude of perihelion to argument of perihelion
      g = p - n
c
c Rotation factors
      call mco_sine (i,si,ci)
      call mco_sine (g,sg,cg)
      call mco_sine (n,sn,cn)
      z1 = cg * cn
      z2 = cg * sn
      z3 = sg * cn
      z4 = sg * sn
      d11 =  z1 - z4*ci
      d12 =  z2 + z3*ci
      d13 = sg * si
      d21 = -z3 - z2*ci
      d22 = -z4 + z1*ci
      d23 = cg * si
c
c Semi-major axis
      a = q / (1.d0 - e)
c
c Ellipse
      if (e.lt.1.d0) then
        romes = sqrt(1.d0 - e*e)
        temp = mco_kep (e,l)
        call mco_sine (temp,se,ce)
        z1 = a * (ce - e)
        z2 = a * romes * se
        temp = sqrt(mu/a) / (1.d0 - e*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
      else
c Parabola
        if (e.eq.1.d0) then
          ce = orbel_zget(l)
          z1 = q * (1.d0 - ce*ce)
          z2 = 2.d0 * q * ce
          z4 = sqrt(2.d0*mu/q) / (1.d0 + ce*ce)
          z3 = -ce * z4
        else
c Hyperbola
          romes = sqrt(e*e - 1.d0)
          temp = orbel_fhybrid(e,l)
          call mco_sinh (temp,se,ce)
          z1 = a * (ce - e)
          z2 = -a * romes * se
          temp = sqrt(mu/abs(a)) / (e*ce - 1.d0)
          z3 = -se * temp
          z4 = romes * ce * temp
        end if
      endif
c
      x = d11*z1 + d21*z2
      y = d12*z1 + d22*z2
      z = d13*z1 + d23*z2
      u = d11*z3 + d21*z4
      v = d12*z3 + d22*z4
      w = d13*z3 + d23*z4
c
c------------------------------------------------------------------------------
c
      return
      end

***********************************************************************
c                    ORBEL_FHYBRID.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           n ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
*	         For larger N, uses FGET
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26,1992.
*     REVISIONS: 
*     REVISIONS: 2/26/93 hfl
***********************************************************************

	real*8 function orbel_fhybrid(e,n)

        implicit NONE    

c...   Version of Swift
       real*8 VER_NUM
       parameter(VER_NUM=2.0d0)

c...   Maximum array size
       integer  NPLMAX, NTPMAX
c       parameter  (NPLMAX = 21)   ! max number of planets, including the Sun
       parameter  (NPLMAX = 51)   ! max number of planets, including the Sun
       parameter  (NTPMAX = 1001) ! max number of test particles

c...   Size of the test particle integer status flag
        integer NSTATP            ! Number of status parameters
        parameter  (NSTATP = 3)
        integer NSTAT            ! Number of status parameters
        parameter  (NSTAT = NSTATP + NPLMAX - 1)  ! include one for @ planet

c...   Size of the test particle integer status flag
        integer NSTATR
        parameter  (NSTATR = NSTAT)  ! io_init_tp assumes NSTAT==NSTATR

c...   convergence criteria for danby
        real*8 DANBYAC , DANBYB
        parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
        integer NLAG1, NLAG2
        parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
        real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
        real*8 PI,TWOPI,PIBY2,DEGRAD
        parameter (PI = 3.14159265358979D0)
        parameter (TWOPI = 2.0D0 * PI)
        parameter (PIBY2 = PI/2.0D0)
        parameter (DEGRAD = 180.0D0 / PI)


c...  Inputs Only: 
	real*8 e,n

c...  Internals:
	real*8 abn
        real*8 orbel_flon,orbel_fget

c----
c...  Executable code 

	abn = n
	if(n.lt.0.d0) abn = -abn

	if(abn .lt. 0.636d0*e -0.6d0) then
	  orbel_fhybrid = orbel_flon(e,n)
	else 
	  orbel_fhybrid = orbel_fget(e,n)
	endif   

	return
	end  ! orbel_fhybrid
c-------------------------------------------------------------------




