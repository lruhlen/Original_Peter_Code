
c..........................................................................

       program nbf


c..........................................................................

c..    N body code: 
c..    Requires Numerical Recipes Routines:  rqks, rk4, ran2,
c...   rkck, rlft3, fourn, gasdev

       include 'fftdefs.h'

       parameter(n1=256,n2=256,n3=1)
       parameter(n=6*nb,pi=3.1415926535)

       dimension y(n),dydx(n),yscal(n)
       dimension iell(ngals)
       dimension sigma(64)

       common /path2/ rm(nb),vcx,vcy,vcz,dorig,eps2,dynam,overscale
       common /particles/ iactive(nb),index(nb,3)
       common /variables/ d(n1,n2,n3),phi(n1,n2,n3)

       external derivs

c..    Fourier cells in each direction
       ngrid=n1
       ngrid2=ngrid/2
       ngridz=1
c..    Overall accuracy criterion
       Accuracy=1.e-2
c..    Timestep accuracy for Integrator
       Acc=1.0e-3
c..    Attempted timestep
       htry=0.05
c..    Softening parameter
       eps2=0.00025
c..    random number initiator
       idum=-2
c..    number of phase space dimensions
       nvar=n
c..    maximum number of timesteps
       nstepsmax=5000
c..    printout interval
       ipint=20
c..    starting time for the integration
       x1=0.
c..    stopping time for the integration
       x2=50.
c..    overall size of the computational box
       overscale=2.0
       o2scale=overscale*2
c..    Initial radius of the galaxy
       galscale=1.0



c.       assign unit mass to galaxy 
         do i=1,nb
           rm(i)=1.0/real(nb)
         end do
         rmtot=1.

         nsnap=0
       
         do i=1,nb

           ib=(i-1)*6

c..        First choose a radius of the current particle.
6108       continue
           do isample=1,nb
             rsamp=galscale*ran2(idum)
             ysamp=ran2(idum)
             if(ysamp.lt.rsamp) then
               goto 6107
             end if
           end do
6107       continue

c..        Next see if current radius makes the surface density cut:
           ysamp=ran2(idum)
           if(ysamp.gt.sqrt(1-rsamp**2)) then
             goto 6108
           end if
c..        Arriving here means particle radius has been determined

c..        random angle:
           fi=2.*pi*ran2(idum)

           y(ib+1)=rsamp*cos(fi)
           y(ib+3)=rsamp*sin(fi)
           y(ib+5)=0.0

c..        set up the circular velocities:
           vphi=0.808*(sqrt(3.*pi)/2.)*rsamp
           y(ib+2)=-vphi*sin(fi)
           y(ib+4)= vphi*cos(fi)
           y(ib+6)= 0.0

c..        set up azimuthal, radial velocity dispersions
           sigphi=0.522*sqrt(1-rsamp**2)
           vsigphi=sigphi*gasdev(idum)

           sigr=sigphi
           vsigr=sigr*gasdev(idum)

           y(ib+2)=y(ib+2)-vsigphi*sin(fi)+vsigr*cos(fi)
           y(ib+4)=y(ib+4)+vsigphi*cos(fi)+vsigr*sin(fi)
           y(ib+6)= 0.0


         end do

c..      check surface density
c         itot=0
c         do i=1,64
c           rzero=real(i-1)/64.
c           rone=real(i)/64.
c           sigma(i)=0.
c           do j=1,nb
c             ib=(j-1)*6
c             if( (sqrt(y(ib+1)**2+y(ib+3)**2).lt.rone).and.
c     +           (sqrt(y(ib+1)**2+y(ib+3)**2).gt.rzero) ) then
c               sigma(i)=sigma(i)+rm(j)
c               itot=itot+1
c             endif
c           end do
c           sigma(i)=sigma(i)/(pi*(rone+rzero)*(rone-rzero))
c           rad=(rone+rzero)/2.
c           write(11,*) rad,sigma(i),(1.5/pi)*sqrt(1-rad**2)
c         end do




c..    compute, subtract off center of mass.
c..    Compute center of mass velocity

       rmtot=0.
       vcx=0.
       vcy=0.
       vcz=0.

       x=x1
       istop=0

       do i=1,nb
         ib=(i-1)*6
         rmtot=rmtot+rm(i)
         vcx=vcx+rm(i)*y(ib+1)
         vcy=vcy+rm(i)*y(ib+3)
         vcz=vcz+rm(i)*y(ib+5)
       end do

       vcx=vcx/rmtot
       vcy=vcy/rmtot
       vcz=vcz/rmtot

       do i=1,nb
         ib=(i-1)*6
         y(ib+1)=y(ib+1)-vcx
         y(ib+3)=y(ib+3)-vcy
         y(ib+5)=y(ib+5)-vcz
       end do

       call printout(x,y,nsnap)
       nsnap=nsnap+1


c..    number of integration steps in current trial
       iflag=0

c..    Compute center of mass velocity

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
    
       do i=1,nb
         iactive(i)=1
       end do

c..    precompute FFT constants:
       call set_gravity3d(overscale,eps2)
       
c-------------------------------------------------------------------------
c..    Start overall integration loop

2105   continue

         iflag=iflag+1
         
c..      Zero mass array
         do i=1,ngrid
           do j=1,ngrid
             do k=1,ngridz
               d(i,j,k)=0.
             end do
           end do
         end do

c..      fill up the mass array
         numb=0
         do ib=1,nb

           if(iactive(ib).eq.1) then

             do j=1,3
               if(j.ne.3) then
                 jb=(ib-1)*6+2*j-1
                 if(abs(y(jb)).gt.overscale) iactive(ib)=0
                 index(ib,j)=int((y(jb)+overscale)*ngrid/o2scale)+1  
               else
                 index(ib,3)=1
               endif
             end do

             numb=numb+1

             d(index(ib,1),index(ib,2),index(ib,3))=
     +       d(index(ib,1),index(ib,2),index(ib,3))+rm(ib)
             if(iactive(ib).eq.1) then

             dtot=dtot+rm(ib)
             end if

           end if

         end do

c..      Compute the potential using the Fourier Method

         call gravity3d

c..      Test calculation of circular velocity:
c         do i=n1/2+1,n1/2+n1/4
c           j=n1/2
c           k=1
c           radius=(real(i-(n1/2))/(n1/4) + real(i-(n1/2)+1)/(n1/4))/2.
c           potential1=(phi((n1/2)+1,j,k)+phi((n1/2)+1,j-1,k))/2.
c           potential=(phi(i,j,k)+phi(i,j-1,k))/2.
c           pot=(3.*pi*radius**2)/8.
c           write(10,*) radius,pot,potential-potential1
c           delx=1./real(n1/4)
c           vc2=radius*(phi(i+1,j,k)-phi(i-1,j,k))/delx
c           vc=sqrt(vc2)
c           write(10,*) i,real(i-(n1/2))/(n1/4),phi(i,j,k),vc,1.53*radius
c            
c         end do

c..      Test calculation of potential:
c         xloc=0.15
c         yloc=0.5
c         zloc=0.
c         potdirect=0.
c         do i=1,nb
c           j1=(i-1)*6+1
c           j2=(i-1)*6+3
c           j3=(i-1)*6+5
c           disti=sqrt((y(j1)-xloc)**2+(y(j2)-yloc)**2+(y(j3)-zloc)**2   
c     +                +eps2)
c           potdirect=potdirect-rm(i)/disti
c         end do
c
c         write(*,*) 'Direct Summation Potential:',potdirect
c
c         ind1=int((xloc+overscale)*ngrid/o2scale)+1
c         ind2=int((yloc+overscale)*ngrid/o2scale)+1
c         ind3=int((zloc+overscale)*ngrid/o2scale)+1
c         ind3=1
c..       fft version of the potential
c         write(*,*) 'FFT Potential:',phi(ind1,ind2,ind3)
c         stop


c..      Use current timestep to take a Runge-Kutta Integration Step
      
         call derivs(x,y,dydx)

         do iscale=1,nvar
           yscal(iscale)=abs(y(iscale))+abs(htry*dydx(iscale))+1.e-30
         end do

c..      Burlirsch-Stoer is unnecessary:
c         call bsstep(y,dydx,nvar,x,htry,acc,yscal,hdid,hnext,derivs)  

c..      Runge Kutta is also a bit of overkill, considering FFT potential:
         call rkqs(y,dydx,nvar,x,htry,acc,yscal,hdid,hnext,derivs)  

          write(*,*) iflag,x,htry,hdid
   
c..      check for printout, 
         if(mod(iflag,ipint).eq.0) then

           call printout(x,y,nsnap)
           nsnap=nsnap+1

c..        check if number of integration steps exceeded:
           if(iflag.gt.nstepsmax) then
             istop=10
             goto 2106
           end if

c..        check if physical time exceeded:
           if(x.gt.x2) then
             istop=11
             goto 2106
           end if

c..        Outcome checklist finished

         end if

         goto 2105

c..      End of overall integration loop
c-------------------------------------------------------------------------

2106     continue

         write(*,*) 'istop',istop
         call printout(x,y,nsnap)

       end




      subroutine derivs(x,y,dydx)

      include 'fftdefs.h'
      parameter(n1=256,n2=256,n3=1)
      parameter(n=6*nb)

      common /path2/ rm(nb),vcx,vcy,vcz,dorig,eps2,dynam,overscale

      common /particles/ iactive(nb),index(nb,3)
      common /variables/ d(n1,n2,n3),phi(n1,n2,n3)


      real*8 x,y(n),dydx(n)

      real*8 denom(ngals,ngals)
      real*8 denom2(nb,ngals)


      do i=2,n,2
        dydx(i-1)=y(i)
      end do

      delx=2.*overscale/real(n1)
      del2x=4.*overscale/real(n1)

      do i=1,nb
        ib=(i-1)*6

        do ic=1,3

          dydx(ib+(2*ic))=0.


c..       fft potential       
          if(iactive(i).eq.1) then
          if (ic.eq.1) then   
            if(index(i,1).eq.1) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
     +        phi(index(i,1)+1,index(i,2),index(i,3))-
     +        phi(index(i,1),index(i,2),index(i,3)))/delx
            elseif(index(i,1).eq.n1) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
     +        phi(index(i,1),index(i,2),index(i,3))-
     +        phi(index(i,1)-1,index(i,2),index(i,3)))/delx
            else
              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
     +        phi(index(i,1)+1,index(i,2),index(i,3))-
     +        phi(index(i,1)-1,index(i,2),index(i,3)))/del2x
            end if
          end if

          if (ic.eq.2) then   
            if(index(i,2).eq.1) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
     +        phi(index(i,1),index(i,2)+1,index(i,3))-
     +        phi(index(i,1),index(i,2),index(i,3)))/delx
            elseif(index(i,2).eq.n2) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
     +        phi(index(i,1),index(i,2),index(i,3))-
     +        phi(index(i,1),index(i,2)-1,index(i,3)))/delx
            else
              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
     +        phi(index(i,1),index(i,2)+1,index(i,3))-
     +        phi(index(i,1),index(i,2)-1,index(i,3)))/del2x
            end if
          end if
         
c          if (ic.eq.3) then   
c            if(index(i,3).eq.1) then
c              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
c     +        phi(index(i,1),index(i,2),index(i,3)+1)-
c     +        phi(index(i,1),index(i,2),index(i,3)))/delx
c            elseif(index(i,3).eq.n3) then
c              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
c     +        phi(index(i,1),index(i,2),index(i,3))-
c     +        phi(index(i,1),index(i,2),index(i,3)-1))/delx
c            else
c              dydx(ib+(2*ic))=dydx(ib+(2*ic))-(
c     +        phi(index(i,1),index(i,2),index(i,3)+1)-
c     +        phi(index(i,1),index(i,2),index(i,3)-1))/del2x
c            end if
c          end if

           end if
  
        end do
      end do


 
      return
      end




c..................................................................
c
      subroutine printout(x,y,nsnap)
c
c..................................................................

      include 'fftdefs.h'


      common /path2/ rm(nb),vcx,vcy,vcz,dorig,eps2,dynam,overscale
      real*8 y(6*nb)

c......local character variables for file manipulation
       character*3 sstring
c       character*29 filestem
       character*4 filestem
       character*10 nstring
c       character*32 target_out
       character*7 target_out


       DATA nstring/'0123456789'/
       save

       sstring(1:1)=nstring(1+nsnap/100:1+nsnap/100)
       istring=1+MOD(nsnap,100)/10
       sstring(2:2)=nstring(istring:istring)
       istring=1+MOD(nsnap,10)
       sstring(3:3)=nstring(istring:istring)

c..    use this if files go to the scratch disk:
c..    filestem= '/scratch/isis/laugh/run3/snap'
c..    target_out= filestem(1:29)//sstring(1:3)

c..    use this if files go to the current directory:
       filestem='snap'
       target_out= filestem(1:4)//sstring(1:3)



       open (unit=10,file=target_out,status='unknown')




c..   write position data to files

      do ib=1,ngals
        jb=(ib-1)*6
        write(10,2104) y(jb+1),y(jb+3),y(jb+5) 
      end do

      do ib=ngals+1,nb
        jb=(ib-1)*6
        write(10,2104) y(jb+1),y(jb+3),y(jb+5) 
      end do

      close(10)

2104       format(4e13.5)
      return
      end




c...................................................................
c
       subroutine set_gravity3d(overscale,eps)
c
c...................................................................

       implicit real*8(a-h,o-z), integer*4(i-n)
       parameter(n1=256,n2=256,n3=1)
       parameter(pi=3.14159)
       common /variables/ d(n1,n2,n3),phi(n1,n2,n3)
       common /potential/ Glm(-n1:n1-1,-n2:n2-1,-n3:n3-1)

       save

c......Prepare fixed terms for calls to FFT 3--D potential solver
c......This routine only called once during run.


       dx=(overscale*2.)/real(n1)
       dy=(overscale*2.)/real(n2)
       dz=(overscale*2.)/real(n3)
       dz=0.0

C..... compute fixed constants for the Poisson solver
       do i=-n1,n1-1
         do j=-n2,n2-1
           do l=-n3,n3-1
             Glm(i,j,l)=-1./sqrt
     +       (eps+(dx*real(i))**2+(dy*real(j))**2+(dz*real(l))**2)  
           end do
         end do
       end do

       return
       end

c...................................................................
c
       subroutine gravity3d
c
c...................................................................

c......Use 3D Fourier Convolution theorem to solve for gravitational
c......potential phi(i,j,l).

       implicit real*8(a-h,o-z), integer*4(i-n)
       parameter(n1=256,n2=256,n3=1)
       parameter(pi=3.14159)
       common /variables/ d(n1,n2,n3),phi(n1,n2,n3)
       common /potential/ Glm(-n1:n1-1,-n2:n2-1,-n3:n3-1)

c......define local variables

       real*8  Mlm(-n1:n1-1,-n2:n2-1,-n3:n3-1)
       real*8 fac,data1(n1*2,n2*2,n3*2),data2(n1*2,n2*2,n3*2)
       complex*16 spec1(n1,n2*2,n3*2),speq1(n2*2,n3*2),
     +            spec2(n1,n2*2,n3*2),speq2(n2*2,n3*2),
     +            zpec1(n1*n2*2*n3*2),zpeq1(n2*2*n3*2),
     +            zpec2(n1*n2*2*n3*2),zpeq2(n2*2*n3*2)

       equivalence (data1,spec1,zpec1), (data2,spec2,zpec2),
     +             (speq1,zpeq1), (speq2,zpeq2)

       save


       do i=-n1,n1-1
         do j=-n2,n2-1
           do l=-n3,n3-1
             if((i.lt.0).or.(j.lt.0).or.(l.lt.0))then
               Mlm(i,j,l)=0.
             else
               Mlm(i,j,l)=d(i+1,j+1,l+1)
             endif
           end do
         end do
       end do

       do i=1,2*n1
         do j=1,2*n2
           do l=1,2*n3
             if     (i.le.n1.and.j.le.n2.and.l.le.n3) then
               data1(i,j,l)=Glm(i-1,j-1,l-1)
               data2(i,j,l)=Mlm(i-1,j-1,l-1)
             elseif (i.gt.n1.and.j.le.n2.and.l.le.n3) then
               data1(i,j,l)=Glm(i-1-2*n1,j-1,l-1)
               data2(i,j,l)=Mlm(i-1-2*n1,j-1,l-1)
             elseif (i.le.n1.and.j.le.n2.and.l.gt.n3) then
               data1(i,j,l)=Glm(i-1,j-1,l-1-2*n3)
               data2(i,j,l)=Mlm(i-1,j-1,l-1-2*n3)
             elseif (i.gt.n1.and.j.le.n2.and.l.gt.n3) then
               data1(i,j,l)=Glm(i-1-2*n1,j-1,l-1-2*n3)
               data2(i,j,l)=Mlm(i-1-2*n1,j-1,l-1-2*n3)
             elseif (i.le.n1.and.j.gt.n2.and.l.gt.n3) then
               data1(i,j,l)=Glm(i-1,j-1-n2*2,l-1-n3*2)
               data2(i,j,l)=Mlm(i-1,j-1-n2*2,l-1-n3*2)
             elseif (i.gt.n1.and.j.gt.n2.and.l.gt.n3) then
               data1(i,j,l)=Glm(i-1-2*n1,j-1-2*n2,l-1-2*n3)
               data2(i,j,l)=Mlm(i-1-2*n1,j-1-2*n2,l-1-2*n3)
             elseif (i.le.n1.and.j.gt.n2.and.l.le.n3) then
               data1(i,j,l)=Glm(i-1,j-1-2*n2,l-1)
               data2(i,j,l)=Mlm(i-1,j-1-2*n2,l-1)
             elseif (i.gt.n1.and.j.gt.n2.and.l.le.n3) then
               data1(i,j,l)=Glm(i-1-2*n1,j-1-2*n2,l-1)
               data2(i,j,l)=Mlm(i-1-2*n1,j-1-2*n2,l-1)
             else
               write(*,*) 'Error in potential grid'
               stop
             end if
           end do
         end do
       end do

c..    call the FFT routines and perform convolution

       call rlft3(data1,speq1,n1*2,n2*2,n3*2,1)
       call rlft3(data2,speq2,n1*2,n2*2,n3*2,1)

       fac=2./(8*N1*N2*N3)
       do j=1,4*n1*n2*n3
         zpec1(j)=fac*zpec1(j)*zpec2(j)
       end do
       do j=1,4*n2*n3
         zpeq1(j)=fac*zpeq1(j)*zpeq2(j)
       end do

       call rlft3(data1,speq1,2*n1,2*n2,2*n3,-1)

       do i=1,n1
         do j=1,n2
           do l=1,n3
             phi(i,j,l)=data1(i,j,l)
           end do
         end do
       end do

       return
       end


       

     





