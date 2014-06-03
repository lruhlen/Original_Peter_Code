

c..     Symplectic map for the 3-body problem 
c..     Reference: Wisdom \& Holman 1991 AJ 102, 1528

      program symplectic

      implicit real*8(a-h,o-z), integer*4(i-n)
      parameter (nstars=1,nb=3)

      parameter(n=6*nb,pi=3.14159265358979323846)

      dimension y(n),yj(n)
      dimension at(nb,nb),ai(nb,nb)
      dimension rm(nb)


       open(unit=1,file='input.symplec',status='unknown')
       open(unit=2,file='planet.2',status='unknown')
       open(unit=3,file='planet.3',status='unknown')

c..    Set physical constants

       rmsun=1.98911e+33
       rau=1.495978707e+13
       G=6.672e-8
       rmjup=1.8986e+30
       year=365.25*24.*3600.

c..    Read in the input data

c..    integration timescale (years)
       read(1,*) x2
c..    timestep size (days)
       read(1,*) tstepsize
c..    Energy accuracy tolerance
       read(1,*) accuracy
c..    printout interval
       read(1,*) iprint


c..    Set the observed parameters (cgs)
       read(1,*) rmstar
       rmstar=rmstar*rmsun
       rm(1)=rmstar

c..    Periods
       read(1,*) period1
       period1=period1*(24.*3600.)
       read(1,*) period2
       period2=period2*(24.*3600.)

       read(1,*) imean
    
       if(imean.eq.0) then

c..      Mean anomoly

         read(1,*) anom1
         read(1,*) anom2

         anom1=((anom1)/360.)*2.*pi
         anom2=((anom2)/360.)*2.*pi

       else

c..      Time of periastron passage:

         read(1,*) Tperi1
         Tperi1=Tperi1*(24.*3600.)
         read(1,*) Tperi2
         Tperi2=Tperi2*(24.*3600.)

       end if


c..    current epoch:
       read(1,*) Epoch
       Epoch=Epoch*(24.*3600.)

c..    If using "Tperi" as input quantities, obtain the mean anomalies
c..    for the two planets:

       if(imean.ne.0) then
         anom1=((2.*pi)/period1)*(epoch-Tperi1)
         anom2=((2.*pi)/period2)*(epoch-Tperi2)
       end if

c..    eccentricities
       read(1,*) ebody1
       read(1,*) ebody2

c..    argument of perihelion
       read(1,*) rOm1
       read(1,*) rOm2

       rOm1=((rOm1)/360.)*2.*pi
       rOm2=((rOm2)/360.)*2.*pi

c..    inclinations
       read(1,*) ri1
       read(1,*) ri2


       ri1=((ri1)/360.)*2.*pi
       ri2=((ri2)/360.)*2.*pi

c..    nodes
       read(1,*) rnode1
       read(1,*) rnode2

       rnode1=((rnode1)/360.)*2.*pi
       rnode2=((rnode2)/360.)*2.*pi
     
       read(1,*) ijov

       if(ijov.eq.1) then
c..      we are inputting Jupiter Masses 

         read(1,*) rm(2)
         rm(2)=rm(2)*rmjup

         read(1,*) rm(3)
         rm(3)=rm(3)*rmjup

         abod1=period1*period1*G*(rmstar+rm(2))
         abod1=abod1/(4.*pi*pi)
         abod1=abod1**0.33333333333

         abod2=period2*period2*G*(rmstar+rm(2)+rm(3))
         abod2=abod2/(4.*pi*pi)
         abod2=abod2**0.33333333333

       else
c..      

         read(1,*) rK1
         rK1=rK1*100.
         read(1,*) rK2
         rK2=rK2*100.

         do i=1,100000000
           rmp=(real(i)/10000.)
           rmp=rmp*rmjup

           abod1=period1*period1*G*(rmstar+rmp)
           abod1=abod1/(4.*pi*pi)
           abod1=abod1**0.33333333333
           q1=(1./rK1)*(2.*pi*G/period1)**0.333333333
           q1=q1/sqrt(1-ebody1*ebody1)
           q1=(q1*rmp*fsini2)/(rmstar+rmp)**0.66666666
           if(q1.gt.1.) then
             rm(2)=rmp
             goto 4104
           end if
         end do
4104     continue

          do i=1,100000000
            rmp=(real(i)/10000.)
            rmp=rmp*rmjup

            abod2=period2*period2*G*(rmstar+rm(2)+rmp)
            abod2=abod2/(4.*pi*pi)
            abod2=abod2**0.33333333333
            q1=(1./rK2)*(2.*pi*G/period2)**0.333333333
            q1=q1/sqrt(1-ebody2*ebody2)
            q1=(q1*rmp*fsini3)/(rmstar+rm(2)+rmp)**0.66666666
            if(q1.gt.1.) then
              rm(3)=rmp
              goto 4105
            end if
          end do
4105      continue
   
        end if

c       write(*,*) 'Semi-Major Axes (in AU) of bodies 1 and 2:'
c       write(*,*) abod1/rau,abod2/rau
c       write(*,*)
c       write(*,*) 'Eccentricities of bodies 1 and 2:'
c       write(*,*)  ebody1,ebody2
c       write(*,*)
c       write(*,*) 'M1=',rm(2)/rmjup,'M2=',rm(3)/rmjup
c       write(*,*)


c..    code uses units G=1,m=1Msun,t=1yr,d=5.091369e+13cm
       rlu=(G*rmsun*(year**2))**0.333333333333

       abod1=abod1/rlu
       abod2=abod2/rlu

       rm(1)=rm(1)/rmsun
       rm(2)=rm(2)/rmsun
       rm(3)=rm(3)/rmsun

       rmf2=rm(1)+rm(2)
       rmf3=rm(1)+rm(3)

       period1=period1/year
       period2=period2/year

c..    set minimum distance flags (to signal close approaches)
       d12min=1.0e+10
       d13min=1.0e+10
       d23min=1.0e+10

c...     Set up cartesian data
         y(1)=0.0
         y(3)=0.0
         y(5)=0.0
         y(2)=0.0
         y(4)=0.0
         y(6)=0.0

         rmu=rm(1)+rm(2)
         q=abod1*(1-ebody1)
         e=ebody1
         rinc=ri1
         p=rOm1
         rn=rnode1
         rl=anom1
         call mco_el2x (rmu,q,e,rinc,p,rn,rl,rx,ry,rz,ru,rv,rw)

         y(7)=rx
         y(9)=ry
         y(11)=rz
         y(8)=ru
         y(10)=rv
         y(12)=rw

         rmu=rm(1)+rm(3)
         q=abod2*(1-ebody2)
         e=ebody2
         rinc=ri2
         p=rOm2
         rn=rnode2
         rl=anom2
         call mco_el2x (rmu,q,e,rinc,p,rn,rl,rx,ry,rz,ru,rv,rw)


         y(13)=rx
         y(15)=ry
         y(17)=rz
         y(14)=ru
         y(16)=rv
         y(18)=rw

         v1=rm(1)+rm(2)
         v2=y(7)-y(1)
         v3=y(9)-y(3)
         v4=y(11)-y(5)
         v5=y(8)-y(2)
         v6=y(10)-y(4)
         v7=y(12)-y(6)
         call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,ri,p,rn,rl)
         q=q/(1-e)

c       write(*,*) 'Inner Planet Starting Conditions:'
c       write(*,*) 'a, e, i (dg), argperi (dg), rn, mean anomaly (dg)'

       ri=(ri/(2.*pi))*360.
       p=  (p/(2.*pi))*360.
       rn=(rn/(2.*pi))*360.
       rl=(rl/(2.*pi))*360.

c       write(*,*) (q*rlu)/rau,e,ri,p,rn,rl

       v1=rm(1)+rm(3)
       v2=y(13)-y(1)
       v3=y(15)-y(3)
       v4=y(17)-y(5)
       v5=y(14)-y(2)
       v6=y(16)-y(4)
       v7=y(18)-y(6)
       call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,ri,p,rn,rl)

       q=q/(1-e)

c       write(*,*) 'Outer Planet Starting Conditions:'
c       write(*,*) 'a, e, i (dg), argperi (dg), rn, mean anomaly (dg)'

       ri=(ri/(2.*pi))*360.
       p=  (p/(2.*pi))*360.
       rn=(rn/(2.*pi))*360.
       rl=(rl/(2.*pi))*360.

c       write(*,*) (q*rlu)/rau,e,ri,p,rn,rl



c..    Determine System Energy

         call EnergySum(rm,y,Energy)

         Eorig=Energy

c..    Various Mass combinations

         rmsum=rm(1)+rm(2)+rm(3)

         rm12=rm(1)+rm(2)

         rm3rm1=rm(3)*rm(1)
         rm3rm2=rm(3)*rm(2)

c..    Eq. 10, WH91
         eta0=rm(1)
         eta1=rm(1)+rm(2)
         eta2=rm(1)+rm(2)+rm(3)

c..    Jacobian mass factors
         rmj0=eta2
         rmj1=eta0*rm(2)/eta1
         rmj2=eta1*rm(3)/eta2

c..    Mass factors sent to Keplerian orbit integrator
         rjkmf1=rm(2)*rm(1)/rmj1
         rjkmf2=rm(3)*rm(1)/rmj2

c..    Mass factors used in Jacobian velocity updates
         rm31d1=rm3rm1/rmj1
         rm32d1=rm3rm2/rmj1
         rm31d2=rm3rm1/rmj2
         rm32d2=rm3rm2/rmj2
 
c..   Transform initial model to Jacobi coordinates

c..   Build transformation matrix for y(n) -> yj(n)

      at(1,1)=rm(1)/rmsum
      at(1,2)=rm(2)/rmsum
      at(1,3)=rm(3)/rmsum
      at(2,1)=-1.
      at(2,2)=1.
      at(2,3)=0.
      at(3,1)=-rm(1)/rm12
      at(3,2)=-rm(2)/rm12
      at(3,3)=1.

c..   Build inverse transformation matrix for yj(n) -> y(n)

      ai(1,1)=1.0
      ai(1,2)=-1.0*rm(2)/rm12 
      ai(1,3)=-1.0*rm(3)/rmsum 
      ai(2,1)=1.0
      ai(2,2)=rm(1)/rm12 
      ai(2,3)=-1.0*rm(3)/rmsum 
      ai(3,1)=1.0
      ai(3,2)=0.0
      ai(3,3)=rm12/rmsum 

c..   Load Jacobi coordinate yj vector:

        yj( 1)=at(1,1)*y( 1)+at(1,2)*y( 7)+at(1,3)*y(13)
        yj( 7)=at(2,1)*y( 1)+at(2,2)*y( 7)+at(2,3)*y(13)
        yj(13)=at(3,1)*y( 1)+at(3,2)*y( 7)+at(3,3)*y(13)

        yj( 3)=at(1,1)*y( 3)+at(1,2)*y( 9)+at(1,3)*y(15)
        yj( 9)=at(2,1)*y( 3)+at(2,2)*y( 9)+at(2,3)*y(15)
        yj(15)=at(3,1)*y( 3)+at(3,2)*y( 9)+at(3,3)*y(15)

        yj( 5)=at(1,1)*y( 5)+at(1,2)*y(11)+at(1,3)*y(17)
        yj(11)=at(2,1)*y( 5)+at(2,2)*y(11)+at(2,3)*y(17)
        yj(17)=at(3,1)*y( 5)+at(3,2)*y(11)+at(3,3)*y(17)

        yj( 2)=at(1,1)*y( 2)+at(1,2)*y( 8)+at(1,3)*y(14)
        yj( 8)=at(2,1)*y( 2)+at(2,2)*y( 8)+at(2,3)*y(14)
        yj(14)=at(3,1)*y( 2)+at(3,2)*y( 8)+at(3,3)*y(14)

        yj( 4)=at(1,1)*y( 4)+at(1,2)*y(10)+at(1,3)*y(16)
        yj(10)=at(2,1)*y( 4)+at(2,2)*y(10)+at(2,3)*y(16)
        yj(16)=at(3,1)*y( 4)+at(3,2)*y(10)+at(3,3)*y(16)

        yj( 6)=at(1,1)*y( 6)+at(1,2)*y(12)+at(1,3)*y(18)
        yj(12)=at(2,1)*y( 6)+at(2,2)*y(12)+at(2,3)*y(18)
        yj(18)=at(3,1)*y( 6)+at(3,2)*y(12)+at(3,3)*y(18)


c..   Set timestep in days:

      delt=tstepsize/365.

c..   initialize the time:
      time=0.0

       i=1

c..    Start time integration loop 
2104   continue

c..    leapfrog time intervals
        if (i.eq.1) then
          dt=delt/2.
        else
          dt=delt
        end if
        i=i+1

c..   Kepler advance planet #1
        ib=2
        call kepler(yj,rjkmf1,dt,ib)

c..   Kepler advance planet #2
        ib=3
        call kepler(yj,rjkmf2,dt,ib)


c===================================================
c	interaction hamiltonian terms for 3 bodies 
c===================================================

c..  Term 1 -- interaction r12: 
	apar=-1.0*rm(1)/rm12 
	coeff=rm(2)*rm(3) 
	xx=yj(13) + apar*yj(7) 
	yy=yj(15) + apar*yj(9) 
	zz=yj(17) + apar*yj(11) 
	rsquare=xx*xx + yy*yy + zz*zz 
	rcube=rsquare*Dsqrt(rsquare) 
	dhdq2a=coeff*apar*xx/rcube 
	dhdq2b=coeff*apar*yy/rcube 
	dhdq2c=coeff*apar*zz/rcube 
	dhdq3x=coeff*xx/rcube 
	dhdq3y=coeff*yy/rcube 
	dhdq3z=coeff*zz/rcube 

c..  Term 2 -- the r20 piece 
	apar=rm(2)/rm12
	coeff=rm(1)*rm(3) 
	xx=yj(13) + apar*yj(7) 
	yy=yj(15) + apar*yj(9) 
	zz=yj(17) + apar*yj(11) 
	rsquare=xx*xx + yy*yy + zz*zz 
	rcube=rsquare*Dsqrt(rsquare) 
	dhdq2a=dhdq2a + coeff*apar*xx/rcube 
	dhdq2b=dhdq2b + coeff*apar*yy/rcube 
	dhdq2c=dhdq2c + coeff*apar*zz/rcube 
	dhdq3x=dhdq3x + coeff*xx/rcube 
	dhdq3y=dhdq3y + coeff*yy/rcube 
	dhdq3z=dhdq3z + coeff*zz/rcube 

c..  Term 3 -- the r2prime piece 
	coeff=rm(1)*rm(3) 
	xx=yj(13) 
	yy=yj(15) 
	zz=yj(17) 
	rsquare=xx*xx + yy*yy + zz*zz 
	rcube=rsquare*Dsqrt(rsquare) 
	dhdq3x=dhdq3x - coeff*xx/rcube 
	dhdq3y=dhdq3y - coeff*yy/rcube 
	dhdq3z=dhdq3z - coeff*zz/rcube 

c..  Update the momenta with the kicks: 

	yj(8)=yj(8) - dt*dhdq2a/rmj1
	yj(10)=yj(10) - dt*dhdq2b/rmj1
	yj(12)=yj(12) - dt*dhdq2c/rmj1
	yj(14)=yj(14) - dt*dhdq3x/rmj2
	yj(16)=yj(16) - dt*dhdq3y/rmj2
	yj(18)=yj(18) - dt*dhdq3z/rmj2

c..     Re-load the ordinary coordinates in order to make proximity test

        y( 1)=ai(1,1)*yj( 1)+ai(1,2)*yj( 7)+ai(1,3)*yj(13)
        y( 7)=ai(2,1)*yj( 1)+ai(2,2)*yj( 7)+ai(2,3)*yj(13)
        y(13)=ai(3,1)*yj( 1)+ai(3,2)*yj( 7)+ai(3,3)*yj(13)

        y( 3)=ai(1,1)*yj( 3)+ai(1,2)*yj( 9)+ai(1,3)*yj(15)
        y( 9)=ai(2,1)*yj( 3)+ai(2,2)*yj( 9)+ai(2,3)*yj(15)
        y(15)=ai(3,1)*yj( 3)+ai(3,2)*yj( 9)+ai(3,3)*yj(15)

        y( 5)=ai(1,1)*yj( 5)+ai(1,2)*yj(11)+ai(1,3)*yj(17)
        y(11)=ai(2,1)*yj( 5)+ai(2,2)*yj(11)+ai(2,3)*yj(17)
        y(17)=ai(3,1)*yj( 5)+ai(3,2)*yj(11)+ai(3,3)*yj(17)

        rminus1=y(7)-y(1)
        rminus2=y(9)-y(3)
        rminus3=y(11)-y(5)
        distance12=rminus1*rminus1+rminus2*rminus2+rminus3*rminus3
        if(distance12.lt.d12min) d12min=distance12

        rminus1=y(13)-y(1)
        rminus2=y(15)-y(3)
        rminus3=y(17)-y(5)
        distance13=rminus1*rminus1+rminus2*rminus2+rminus3*rminus3
        if(distance13.lt.d13min) d13min=distance13

        rminus1=y(13)-y(7)
        rminus2=y(15)-y(9)
        rminus3=y(17)-y(11)
        distance23=rminus1*rminus1+rminus2*rminus2+rminus3*rminus3
        if(distance23.lt.d23min) d23min=distance23


        If(mod(i,iprint).eq.0) then


c..       Reload the ordinary velocities:
          y( 2)=ai(1,1)*yj( 2)+ai(1,2)*yj( 8)+ai(1,3)*yj(14)
          y( 8)=ai(2,1)*yj( 2)+ai(2,2)*yj( 8)+ai(2,3)*yj(14)
          y(14)=ai(3,1)*yj( 2)+ai(3,2)*yj( 8)+ai(3,3)*yj(14)

          y( 4)=ai(1,1)*yj( 4)+ai(1,2)*yj(10)+ai(1,3)*yj(16)
          y(10)=ai(2,1)*yj( 4)+ai(2,2)*yj(10)+ai(2,3)*yj(16)
          y(16)=ai(3,1)*yj( 4)+ai(3,2)*yj(10)+ai(3,3)*yj(16)

          y( 6)=ai(1,1)*yj( 6)+ai(1,2)*yj(12)+ai(1,3)*yj(18)
          y(12)=ai(2,1)*yj( 6)+ai(2,2)*yj(12)+ai(2,3)*yj(18)
          y(18)=ai(3,1)*yj( 6)+ai(3,2)*yj(12)+ai(3,3)*yj(18)

c..       Check conservation of energy 
          call EnergySum(rm,y,Energy)
          Efrac=Energy/Eorig
          test=1.0 - Efrac 
          check=Dabs(test)


c..       evaluate orbital elements:
          v1=rm(1)+rm(2)
          v2=y(7)-y(1)
          v3=y(9)-y(3)
          v4=y(11)-y(5)
          v5=y(8)-y(2)
          v6=y(10)-y(4)
          v7=y(12)-y(6)
          call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,ri,p,rn,rl)
          et1=e
          q=q/(1-e)
          ric=ri
          rnc=rn
          rpo=p
      
          q=(q*rlu)/rau

          write(2,2234) time,q*(1-e),q*(1+e),e,ri,p,rn,rl,
     +                d12min,d23min, check
          write(*,2234) time,q*(1-e),q*(1+e),e,ri,p,rn,rl,
     +                d12min,d23min,check
          v1=rm(1)+rm(3)
          v2=y(13)-y(1)
          v3=y(15)-y(3)
          v4=y(17)-y(5)
          v5=y(14)-y(2)
          v6=y(16)-y(4)
          v7=y(18)-y(6)
          call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,ri,p,rn,rl)
          q=q/(1-e)
          q=(q*rlu)/rau

          phicd=dcos(ric)*dcos(ri)+dsin(ric)*dsin(ri)*
     +       cos(rn-rnc)
          phicd=acos(phicd)
          apsided=360.*((p-rpo)/(2.*pi))

          if (apsided.gt.180.) apsided=apsided-360.
          if (apsided.lt.-180.) apsided=apsided+360.


          write(3,2234) time,q*(1-e),q*(1+e),e,ri,p,rn,rl,
     +                apsided,(phicd/(2.*pi))*360.


 2234    format(1x,1Pe14.6,10e12.4)

         end if
c..      end of periodic printout

         time=time+delt
 
         if(time.gt.x2) then
           goto 2106
         end if


       goto 2104
c..    end of timestep loop

2106   continue

         
      end
 

      subroutine kepler(y,rmu,dt,ib)

      implicit real*8(a-h,o-z), integer*4(i-n)
      parameter (nstars=1,nb=3)


      parameter(n=6*nb,pi=3.141592653589793)

      dimension y(n)

      jb=(ib-1)*6
      jb1=jb+1
      jb2=jb+2
      jb3=jb+3
      jb4=jb+4
      jb5=jb+5
      jb6=jb+6

      r0=sqrt(y(jb1)*y(jb1)+y(jb3)*y(jb3)+y(jb5)*y(jb5))
      v0s=y(jb2)*y(jb2)+y(jb4)*y(jb4)+y(jb6)*y(jb6)
      u=y(jb1)*y(jb2)+y(jb3)*y(jb4)+y(jb5)*y(jb6)
      a=1./(2./r0-v0s/rmu)
      en=sqrt(rmu/(a*a*a))
      ec=1.-r0/a
      es=u/(en*a*a)
      e=sqrt(ec*ec + es*es)

c..   solve kepler's equation:
c..   dt passed to this subroutine must be less than one orbit, to avoid:

      dt=dt-nint(en*dt/(2.*pi))*(2.*pi)/en

      yval=en*dt-es
      test=es*cos(yval)+ec*sin(yval)
      if(test.gt.0.) then
        isigma=1
      else
        isigma=-1
      end if
      x=yval+isigma*0.85*e
      nc=0
      del=1.0e-15

c..   start iterative newton's method loop:
2102  continue
        s=sin(x)
        c=cos(x)
        ecs=ec*s
        ecc=ec*c
        ess=es*s
        esc=es*c
        f=x-ecs+es-esc-en*dt
        fp=1.-ecc+ess
        fpp=ecs+esc
        fppp=ecc-ess
        dx=-f/fp
        dx=-f/(fp+dx*fpp/2.)
        dx=-f/(fp+dx*fpp/2.+dx*dx*fppp/6.)
        x=x+dx
        if(abs(dx).lt.del) goto 2104
        nc=nc+1
        if(nc.gt.12) goto 2103
      goto 2102
2103  continue
      write(*,*) 'maximum iteration count in kepler'
2104  continue
      f=(a/r0)*(c-1.)+1.
      g=dt+(s-x)/en
      fdot=-(a/(r0*fp))*en*s
      gdot=(c-1.)/fp + 1.

      yjb1=y(jb1)
      yjb3=y(jb3)
      yjb5=y(jb5)

      y(jb1) = y(jb1)*f + y(jb2)*g
      y(jb3) = y(jb3)*f + y(jb4)*g
      y(jb5) = y(jb5)*f + y(jb6)*g

      y(jb2) = yjb1*fdot +y(jb2)*gdot
      y(jb4) = yjb3*fdot +y(jb4)*gdot
      y(jb6) = yjb5*fdot +y(jb6)*gdot

      return
      end


c...................................................................

      subroutine EnergySum(rm,y,energy)

c...................................................................

      implicit real*8(a-h,o-z), integer*4(i-n)
      parameter (nstars=1,nb=3)


      real*8 y(6*nb),rm(nb)

      dimension denom(nb,nb)
      do i=1,nb
        do j=1,nb
          denom(i,j)=1.e+20
        end do
      end do
      do i=1,nb
        ib=(i-1)*6
        do j=1,nb
          jb=(j-1)*6
          if (i.ne.j) then
            denom(i,j) = (y(jb+1)-y(ib+1))*(y(jb+1)-y(ib+1))
     +                  +(y(jb+3)-y(ib+3))*(y(jb+3)-y(ib+3))
     +                  +(y(jb+5)-y(ib+5))*(y(jb+5)-y(ib+5))
            denom(i,j)=sqrt(denom(i,j))
          end if
        end do
      end do

      Energy=0.
      do i=1,nb
        ib=(i-1)*6
        do j=1,nb
          if(i.ne.j) then
            Energy=Energy-0.5*rm(i)*rm(j)/denom(i,j)
          end if
        end do
        do ic=1,3
          Energy=Energy+0.5*rm(i)*y(ib+2*ic)*y(ib+2*ic)
        end do
      end do

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



