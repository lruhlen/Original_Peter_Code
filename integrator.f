
c..    n-planet integrator, with restart capability

c..    Requires Numerical Recipes Routines: bsstep, ran2, pzextr, mmid


       program bulst

       implicit real*8(a-h,o-z), integer*4(i-n)

       parameter (nb=10)
       parameter(n=6*nb,pi=3.14159265358979323846)


       dimension y(n),dydx(n),yscal(n)
       dimension abod(nb),ebody(nb),ri(nb),rOm(nb),fsini(nb)
       dimension period(nb),rnode(nb),anom(nb),Tperi(nb),rK(nb)
    

       common /path2/ rm(nb)

       external derivs


c..    Set physical constants

       rmsun=1.98911e+33
       rau=1.495978707e+13
       G=6.672e-8
       rmjup=1.8986e+30
       year=365.25*24.*3600.


c..    open the input file
       open(unit=1,file='input',status='old')

c..    open the output files (one for each planet)
c..    (only nb-1 files will be used)

       open(unit=11,file='planet.1',status='unknown')
       open(unit=12,file='planet.2',status='unknown')
       open(unit=13,file='planet.3',status='unknown')
       open(unit=14,file='planet.4',status='unknown')
       open(unit=15,file='planet.5',status='unknown')
       open(unit=16,file='planet.6',status='unknown')
       open(unit=17,file='planet.7',status='unknown')
       open(unit=18,file='planet.8',status='unknown')
       open(unit=19,file='planet.9',status='unknown')

c..    read in start-up conditions from the input file
c..    integration duration (years)
       read(1,*) x2
c..    printout interval
       read(1,*) nprint
c..    mass of the central star
       read(1,*) rmstar
       rmstar=rmstar*rmsun
       rm(1)=rmstar

c..    do we produce a surface of section?
       read(1,*) isection
       if(isection.eq.1) then
         open(unit=20,file='section.data',status='unknown')
         icross=0
         ncross=0
       end if

       read(1,*) iperiod

c..    Periods or semi-major axes
       do i=2,nb

         if (iperiod.eq.0) then
           read(1,*) abod(i)
           abod(i)=abod(i)*rau
         else
           read(1,*) period(i)
           period(i)=period(i)*(24.*3600.)
         end if
         
       end do

c..    Mean anomaly (0), T Peri (1), or mean longitude (2)
       read(1,*) imean
       if(imean.eq.0) then
         do i=2,nb
           read(1,*) anom(i)
           anom(i)=((anom(i))/360.)*2.*pi
         end do
       elseif(imean.eq.1) then
c..      alternately read time of periastron passage
         do i=2,nb
           read(1,*) Tperi(i)
           Tperi(i)=Tperi(i)*(24.*3600.)
         end do
       elseif(imean.eq.2) then
         do i=2,nb
           read(1,*) anom(i)
           anom(i)=((anom(i))/360.)*2.*pi
         end do
       end if

c..    current epoch:
       read(1,*) Epoch
       epochdays=epoch
       Epoch=Epoch*(24.*3600.)

c..    Observational value of radial velocity at first observational epoch.
c..    This will only be used if running comparison curves through RV data
       read(1,*) radVel1

c..    eccentricities
       do i=2,nb
         read(1,*) ebody(i)
       end do

c..    argument of perihelion
       do i=2,nb
         read(1,*) rOm(i)
         rOm(i)=((rOm(i))/360.)*2.*pi
       end do

c..    adjust time of perihelion passage or mean longitude to mean anomaly
       if(imean.eq.1) then
         do i=2,nb
           anom(i)=((2.*pi)/period(i))*(epoch-Tperi(i))
         end do
       end if
       if(imean.eq.2) then
         do i=2,nb
           anom(i)=anom(i)-rOm(i)
         end do
       end if

c..    inclinations
       do i=2,nb
         read(1,*) ri(i)
         fsini(i)=sin(((90.-abs(ri(i)))/360.)*2.*pi)
         ri(i)=((ri(i))/360.)*2.*pi
       end do

c..    nodes
       do i=2,nb
         read(1,*) rnode(i)
         rnode(i)=((rnode(i))/360.)*2.*pi
       end do


c..    if ijov=1, then we are reading in masses directly
       read(1,*) ijov

       if(ijov.eq.1) then

         do i=2,nb
           read(1,*) rm(i)
           rm(i)=rm(i)*(1.0e+27)
           rmsum=0.
           do j=1,i
             rmsum=rmsum+rm(j)
           end do
           if(iperiod.eq.1) then
             abod(i)=period(i)*period(i)*G*(rmsum)
             abod(i)=abod(i)/(4.*pi*pi)
             abod(i)=abod(i)**0.33333333333
           end if
           if(iperiod.eq.0) then
             period(i)=sqrt(abod(i)**3*4*pi*pi/(G*rmsum))
           end if
         end do

c..    if ijov=0, then we get the masses from radial velocity half-amplitudes
       elseif(ijov.eq.0) then

         rminterior=0.
         do j=2,nb
           read(1,*) rK(j)
           rK(j)=rK(j)*100.

           do i=1,100000000
             rmp=(real(i)/10000.)
             rmp=rmp*rmjup

             abod(j)=period(j)*period(j)*G*(rmstar+rmp+rminterior)
             abod(j)=abod(j)/(4.*pi*pi)
             abod(j)=abod(j)**0.33333333333
             q1=(1./rK(j))*(2.*pi*G/period(j))**0.333333333
             q1=q1/sqrt(1-ebody(j)*ebody(j))
             q1=(q1*rmp*fsini(j))/(rmstar+rmp+rminterior)**0.66666666
             if(q1.gt.1.) then
               rm(j)=rmp
               goto 4104
             end if
           end do
4104       continue
           rminterior=rminterior+rm(j)
         end do
       end if


c..    code uses units G=1,m=1Msun,t=1yr,d=5.091369e+13cm
       rlu=(G*rmsun*(year**2))**0.333333333333

       do i=2,nb
         abod(i)=abod(i)/rlu
       end do

       do i=1,nb
         rm(i)=rm(i)/rmsun
       end do

       do i=2,nb
         period(i)=period(i)/year
       end do

c..    Timestep accuracy for Bulirsch-Stoer
       read(1,*) acc

c..    timestep length for integrator (fraction of Period(2))
       read(1,*) hfactor
       hstep=hfactor*period(2)

c..    Jacobi (ijac=1) or astrometric (ijac=0) coordinates
       read(1,*) ijac


c..    miscellaneous initializations
c..    number of phase space dimensions
       nvar=n
c..    starting time for the integration
       x=0.
c..    number of integration steps completed
       iflag=0
         
         
c..    convert initial orbital elements to cartesian initial conditions in star-centered frame:
c..    options for astrocentric or jacobi coordinates:

       if(ijac.eq.0) then

         y(1)=0.0
         y(3)=0.0 
         y(5)=0.0 
         y(2)=0.0
         y(4)=0.0 
         y(6)=0.0 

         do i=2,nb

           rmu=rm(1)+rm(i)
           q=abod(i)*(1.-ebody(i))
           e=ebody(i)
           rinc=ri(i)
           p=rOm(i)
           rn=rnode(i)
           rl=anom(i)

c..        this routine does the elements to cartesian conversion
           call mco_el2x (rmu,q,e,rinc,p,rn,rl,rx,ry,rz,ru,rv,rw)

           y(6*(i-1)+1)=rx
           y(6*(i-1)+3)=ry
           y(6*(i-1)+5)=rz
           y(6*(i-1)+2)=ru
           y(6*(i-1)+4)=rv
           y(6*(i-1)+6)=rw

         end do
   
       elseif (ijac.eq.1) then


         y(1)=0.0
         y(3)=0.0 
         y(5)=0.0 
         y(2)=0.0
         y(4)=0.0 
         y(6)=0.0 

         do i=2,nb

           xcom=0.
           ycom=0.
           zcom=0.
           vxcom=0.
           vycom=0.
           vzcom=0.

           do jdum=1,i-1 
             xcom=xcom+rm(jdum)*y(6*(jdum-1)+1)
             ycom=ycom+rm(jdum)*y(6*(jdum-1)+3)
             zcom=zcom+rm(jdum)*y(6*(jdum-1)+5)
             vxcom=vxcom+rm(jdum)*y(6*(jdum-1)+2)
             vycom=vycom+rm(jdum)*y(6*(jdum-1)+4)
             vzcom=vzcom+rm(jdum)*y(6*(jdum-1)+6)
           end do

           rmu=0.0
           do jdum=1,i-1
             rmu=rmu+rm(jdum)
           end do

           xcom=xcom/rmu
           ycom=ycom/rmu
           zcom=zcom/rmu
           vxcom=vxcom/rmu
           vycom=vycom/rmu
           vzcom=vzcom/rmu

           rmu=rmu+rm(i)

           q=abod(i)*(1.-ebody(i))
           e=ebody(i)
           rinc=ri(i)
           p=rOm(i)
           rn=rnode(i)
           rl=anom(i)

c..        this routine does the elements to cartesian conversion
           call mco_el2x (rmu,q,e,rinc,p,rn,rl,rx,ry,rz,ru,rv,rw)

           y(6*(i-1)+1)=rx+xcom
           y(6*(i-1)+3)=ry+ycom
           y(6*(i-1)+5)=rz+zcom
           y(6*(i-1)+2)=ru+vxcom
           y(6*(i-1)+4)=rv+vycom
           y(6*(i-1)+6)=rw+vzcom

         end do
      
       end if
   

c..    evaluate the orbital elements
       if(ijac.eq.0)then

         do i=2,nb
           v1=rm(1)+rm(i)
           v2=y(6*(i-1)+1)-y(1)
           v3=y(6*(i-1)+3)-y(3)
           v4=y(6*(i-1)+5)-y(5)
           v5=y(6*(i-1)+2)-y(2)
           v6=y(6*(i-1)+4)-y(4)
           v7=y(6*(i-1)+6)-y(6)

           call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,rinc,p,rn,rl)  

           q=q/(1-e)
           inb1 = i-1
           write(*,*) 'Planet ',inb1 , '  starting Conditions:'
           write(*,*) 'a, e, i (dg), argper (dg), rn, m. anom. (dg)'

           rinc=(rinc/(2.*pi))*360.
           p=  (p/(2.*pi))*360.
           rn=(rn/(2.*pi))*360.
           rl=(rl/(2.*pi))*360.

           write(*,*) (q*rlu)/rau,e,rinc,p,rn,rl

         end do

       elseif (ijac.eq.1) then

         do i=2,nb
           xcom=0.
           ycom=0.
           zcom=0.
           vxcom=0.
           vycom=0.
           vzcom=0.

           do jdum=1,i-1 
             xcom=xcom+rm(jdum)*y(6*(jdum-1)+1)
             ycom=ycom+rm(jdum)*y(6*(jdum-1)+3)
             zcom=zcom+rm(jdum)*y(6*(jdum-1)+5)
             vxcom=vxcom+rm(jdum)*y(6*(jdum-1)+2)
             vycom=vycom+rm(jdum)*y(6*(jdum-1)+4)
             vzcom=vzcom+rm(jdum)*y(6*(jdum-1)+6)
           end do

           rmu=0.0
           do jdum=1,i-1
             rmu=rmu+rm(jdum)
           end do

           xcom=xcom/rmu
           ycom=ycom/rmu
           zcom=zcom/rmu
           vxcom=vxcom/rmu
           vycom=vycom/rmu
           vzcom=vzcom/rmu

           rmu=rmu+rm(i)

           v1=rmu
           v2=y(6*(i-1)+1)-xcom
           v3=y(6*(i-1)+3)-ycom
           v4=y(6*(i-1)+5)-zcom
           v5=y(6*(i-1)+2)-vxcom
           v6=y(6*(i-1)+4)-vycom
           v7=y(6*(i-1)+6)-vzcom

           call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,rinc,p,rn,rl)  

           q=q/(1-e)

           write(*,*) 'Planet ',i,'  starting Conditions:'
           write(*,*) 'a, e, i (dg), argper (dg), rn, m. anom. (dg)'

           rinc=(rinc/(2.*pi))*360.
           p=  (p/(2.*pi))*360.
           rn=(rn/(2.*pi))*360.
           rl=(rl/(2.*pi))*360.

           write(*,*) (q*rlu)/rau,e,rinc,p,rn,rl

         end do

      end if




c..    Compute and subtract off center-of-mass velocity

       rmtot=0.
       vcx=0.
       vcy=0.
       vcz=0.

       do i=1,nb
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
         
c..    Determine initial System Energy
       call EnergySum(y,Energy)
       Eorig=Energy

c..    Initial radial velocity of the central star
       vel1=-(y(4)*(rlu/year))/100.



c..    Initializations now finished. 
c..    Start the overall integration loop for the system
2105   continue

         iflag=iflag+1

         htry=hstep

c..      Use current timestep to take a Bulirsch-Stoer Integration Step
         call derivs(x,y,dydx)

c..      refer accuracy to expected phase space values:
         do iscale=1,nvar
           yscal(iscale)=abs(y(iscale))+abs(htry*dydx(iscale))+1.e-30
         end do

         call bsstep(y,dydx,nvar,x,htry,acc,yscal,hdid,hnext,derivs)  

     
c..      Probe the system at cadence nprint, or alternatively, check for axis crossing if isection=1:
         if(mod(iflag,nprint).eq.0.or.isection.eq.1) then

c..        Check conservation
           call EnergySum(y,Energy)
           Efrac=Energy/Eorig
           check=dabs(1.-Efrac)

      
c..        evaluate orbital elements using John Chamber's Routine:

           if(mod(iflag,nprint).eq.0) then
             write(*,*) x

             if(ijac.eq.0)then

               do i=2,nb

                 v1=rm(1)+rm(i)
                 v2=y(6*(i-1)+1)-y(1)
                 v3=y(6*(i-1)+3)-y(3)
                 v4=y(6*(i-1)+5)-y(5)
                 v5=y(6*(i-1)+2)-y(2)
                 v6=y(6*(i-1)+4)-y(4)
                 v7=y(6*(i-1)+6)-y(6)

                 call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,rinc,p,rn,rl)  

                 q=q/(1-e)
                 rinc=(rinc/(2.*pi))*360.
                 p=  (p/(2.*pi))*360.
                 rn=(rn/(2.*pi))*360.
                 rl=(rl/(2.*pi))*360.
                 ifile=9+i

                 write(ifile,2134) x,(q*rlu)/rau,e,rinc,p,rn,rl
2134             format(1p7e13.6)

               end do

             elseif (ijac.eq.1) then

               do i=2,nb

                 xcom=0.
                 ycom=0.
                 zcom=0.
                 vxcom=0.
                 vycom=0.
                 vzcom=0.

                 do jdum=1,i-1 
                   xcom=xcom+rm(jdum)*y(6*(jdum-1)+1)
                   ycom=ycom+rm(jdum)*y(6*(jdum-1)+3)
                   zcom=zcom+rm(jdum)*y(6*(jdum-1)+5)
                   vxcom=vxcom+rm(jdum)*y(6*(jdum-1)+2)
                   vycom=vycom+rm(jdum)*y(6*(jdum-1)+4)
                   vzcom=vzcom+rm(jdum)*y(6*(jdum-1)+6)
                 end do

                 rmu=0.0
                 do jdum=1,i-1
                   rmu=rmu+rm(jdum)
                 end do

                 xcom=xcom/rmu
                 ycom=ycom/rmu
                 zcom=zcom/rmu
                 vxcom=vxcom/rmu
                 vycom=vycom/rmu
                 vzcom=vzcom/rmu

                 rmu=rmu+rm(i)

                 v1=rmu
                 v2=y(6*(i-1)+1)-xcom
                 v3=y(6*(i-1)+3)-ycom
                 v4=y(6*(i-1)+5)-zcom
                 v5=y(6*(i-1)+2)-vxcom
                 v6=y(6*(i-1)+4)-vycom
                 v7=y(6*(i-1)+6)-vzcom

                 call mco_x2el(v1,v2,v3,v4,v5,v6,v7,q,e,rinc,p,rn,rl)  

                 q=q/(1-e)
                 rinc=(rinc/(2.*pi))*360.
                 p=  (p/(2.*pi))*360.
                 rn=(rn/(2.*pi))*360.
                 rl=(rl/(2.*pi))*360.

                 ifile=i+9
                 write(ifile,2134) x,(q*rlu)/rau,e,rinc,p,rn,rl

               end do
 
             end if

           end if

c..        check for periastron passage of second planet in order to construct surface of section
c..        This block occurs if 
           if(isection.eq.1) then

             if(rl.lt.0.5 .and. icross.eq.0) then

               icross=1
   
               if((y(13).gt.0.).and.(y(15).gt.0.)) then
                 phase3=atan(y(15)/y(13))
               elseif((y(13).lt.0.).and.(y(15).gt.0.)) then
                 phase3=atan(y(15)/y(13))+pi
               elseif((y(13).lt.0.).and.(y(15).lt.0.)) then
                 phase3=atan(y(15)/y(13))+pi
               else
                 phase3=atan(y(15)/y(13))+2*pi
               end if

               if((y(7).gt.0.).and.(y(9).gt.0.)) then
                 phase2=atan(y(9)/y(7))
               elseif((y(7).lt.0.).and.(y(9).gt.0.)) then
                 phase2=atan(y(9)/y(7))+pi
               elseif((y(7).lt.0.).and.(y(9).lt.0.)) then
                 phase2=atan(y(9)/y(7))+pi
               else
                 phase2=atan(y(9)/y(7))+2*pi
               end if

               diff=phase3-phase2
               ratio=q3/q2

               if(diff.lt.0.) diff=diff+2.*pi
               ncross=ncross+1
  
               if(mod(ncross,1).eq.0) then
                 write(20,*) diff,ratio,x
               end if

             end if

             if(rl.gt.2. .and. icross.eq.1) then
               icross=0
             end if
   
           end if
c..        end surface of section calculation

2234       format(1p8e11.4)


c..          check if physical time exceeded:
             if(x.gt.x2) then

               istop=2
               goto 2106

             end if

           end if

c..      return to beginning of integration loop
         goto 2105

c..      we've finished up. Write the finish code to 'finish'
2106     continue

         open(unit=30,file='finish',status='unknown')
         write(30,*)  istop
         close(30)

c..      and we're done!


       end




      subroutine derivs(x,y,dydx)

      implicit real*8(a-h,o-z), integer*4(i-n)
      parameter (nb=10)


      common /path2/ rm(nb)

      parameter(n=6*nb)

      real*8 x,y(n),dydx(n)
      real*8 denom(nb,nb)

      

      do i=2,n,2
        dydx(i-1)=y(i)
      end do

      do i=1,nb
        do j=1,nb
          denom(i,j)=1.
        end do
      end do


      do i=1,nb
        do j=i+1,nb

          if (i.ne.j) then

            jb=(j-1)*6
            ib=(i-1)*6

            ytx=y(jb+1)-y(ib+1)
            yty=y(jb+3)-y(ib+3)
            ytz=y(jb+5)-y(ib+5)

            denom(i,j)=ytx*ytx + yty*yty + ytz*ytz
            denom(i,j)=sqrt(denom(i,j))*denom(i,j)
            denom(j,i)=denom(i,j)

          end if

        end do
      end do

      do i=1,nb
        ib=(i-1)*6
        do ic=1,3
          dydx(ib+(2*ic))=0.
          do j=1,nb
            jb=(j-1)*6
            if( i.ne.j) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic)) -
     +        rm(j)*(y(ib+2*ic-1)-y(jb+2*ic-1))/denom(i,j)
            end if
          end do
        end do
      end do


      return
      end


c...................................................................

      subroutine EnergySum(y,energy)

c...................................................................

      implicit real*8(a-h,o-z), integer*4(i-n)
      parameter (nb=10)


      common /path2/ rm(nb)
      real*8 y(6*nb)

      dimension denom(nb,nb)
      do i=1,nb
        do j=1,nb
          denom(i,j)=1.e+10
        end do
      end do
      do i=1,nb
        ib=(i-1)*6
        do j=1,nb
          jb=(j-1)*6
          if (i.ne.j) then
            denom(i,j) = (y(jb+1)-y(ib+1))**2
     +                 +(y(jb+3)-y(ib+3))**2
     +                 +(y(jb+5)-y(ib+5))**2
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
          Energy=Energy+0.5*rm(i)*y(ib+2*ic)**2
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




