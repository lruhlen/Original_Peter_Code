       subroutine nucrat(Rho,T,np,nd,nhe3,dtime,E,                      &
     &                  Rpp,Rhe3he3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,      &
     &                  Xhe4,Xc,Xn,Xo,xpp23,iflip,iwrite)
      include 'parm.h'
c......
       double precision k,np,nde,nd,N0,kb,Lambda0,mu,nhe3
       double precision n(14),nc,nn,no
       COMMON/ABUND/XCA(14),H1(14),A(14)

       dimension Xz(14),AA(14),XBA(14)
       data Xz/1.d0,2.d0,6.d0,11.d0,12.d0,13.d0,14.d0,16.d0,19.d0,20.d0,&
     &        26.d0,7.d0,8.d0,10.d0/
       DATA AA/1.0081451d0,4.003874d0,12.0038156d0,23.d0,24.32d0,       &
     &        26.97d0,28.06d0,32.07d0,39.102d0,40.08d0,55.85d0,         &
     &        14.0075257d0,16.d0,19.99d0/
       DATA Z0,Z1,Z13,HMASS/0.d0,1.d0,0.333333333333333333d0,1.6732D-24/
       DATA CV,CV2,CA,CA2/.934d0,.8724d0,0.5d0,0.25d0/
       DATA ifirst/0/

c......recalculate helium abundance


c......Given the temperature, the density and an abundance set,
c......this routine returns the nuclear energy generation rate
c......in [ergs/sec/gm]. The reactions considered are:

c......proton:proton                  1H(p,e+ve)2H
c......deuterium:proton               2H(p,y)3He
c......He3:He3                        3He(3He,H+H)4He
c......Corrections for ppII & ppIII
c......Equilibrium CNO
c......He4:He4:He4                    4He(aa,y)12C
c......C:He4                          12C(a,y)16O
c......Various high T cooling processes

       save
       j = abs(iwrite)
       if(ifirst.eq.0) then
         ifirst=1
         do i=1,14
           A(i)=AA(i)
           xba(i)=xca(i)
         enddo
c        xmetal=x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(14)
       endif
c
       Rpp = Z0
       Xpd = Z0
       Rhe3he3=Z0
       Rcno=Z0
       Xpc=Z0
       Xpn=Z0
       Xpo=Z0
       R3a=Z0
       R12a=Z0
       xpp23=Z1
       E=Z0
       if(T.lt.5.d5) return
c
c......recompute mass fractions
       xh=(np/rho)*HMASS
       xhe3=3.*(nhe3/rho)*HMASS
       xba(1) = xh/a(1) 
       xba(2) = xhe4/a(2) 
       xba(3) = xc/a(3) 
       xba(12)= xn/a(12) 
       xba(13)= xo/a(13) 
       fpd=Z1
       Zd=Z1
       Zp=Z1
       Zhe3 = 2.d0
       k=1.380658e-16

       zeta=Z0
       do i=1,14
         zeta = zeta+(Xz(i)**2+Xz(i))*xba(i)
       end do
C*******  
       if(T.ge.3.d7) then
          T8=T*1.d-8
C  f3a formula from Clayton (1968)
          f3a=exp(2.76d-3*sqrt(RHO/T8)/T8)
C  Triple alpha from Kippenhahn & Weigert (1991)
          E3a=5.1d11 * f3a * RHO**2 * (Xhe4/T8)**3 * exp(-44.027d0/T8)
          Q3a=1.164d-5
          R3a=E3a*RHO/Q3a
C  Carbon:alpha from Kippenhahn & Weigert(1991)
          T7=T8*10.d0
          f12a=exp(0.071d0*sqrt(zeta*RHO/T7**3))
          T823=T8**(2.d0/3.d0)
          g12a=((Z1+0.134d0*T823)/(Z1+0.01d0*T823))**2
          E12a=1.3d27 * f12a * Xc * Xhe4 * RHO / T8**2 * g12a *
     *         exp(-69.20d0/T8**Z13)
          Q12a=1.146d-5
          R12a=E12a * RHO / Q12a
          E  =E3a + E12a
C Neutrino loss rates Beaudet, Petrosian, and Salpeter  1967
          T9=T8/10.d0
          alam=T*1.683d-10
          axi=((rho/2.d0/1.d9)**(.33333333d0))/alam
          fpair1 = (6.002d19 + 2.084d20*axi + 1.872d21*axi**2)
     &     *exp(-5.5924*axi)
          fpair2=axi**3 + 9.383d-1/alam - 4.141d-1/alam/alam
     &     + 5.829d-2/alam**3
          fpair=fpair1/fpair2
          alam2=alam*alam
          alam4=alam2*alam2
          alam6=alam4*alam2
          alam8=alam6*alam2
          glam=1.d0 -13.04d0*alam2+133.5d0*alam4 
     &    + 1534.d0*alam6 + 918.6d0*alam8 
          epair=glam/rho*fpair*exp(-2.d0/alam)
C         Z2A=xhe4 + xc*3.d0 + xo*4.d0
C         ebrems=0.76d0*Z2A*T8**6
C         if(rho.lt. 1.d8) then
C         factor=2.5d0*(8.d0- log10(rho)) -1.d0
C    &      +0.25d0*log10(rho)
C         ebrems = ebrems/factor
C         end if 
          ebrems=0.
          fphot1 = (4.886d10 + 7.580d10*axi + 6.023d10*axi**2)
     &     *exp(-1.5654*axi)
          fphot2=axi**3 + 6.29d-3/alam + 7.483d-3/alam/alam
     &     + 3.061d-4/alam**3
          fphot=fphot1/fphot2
C    ********  mue=2  *****************
          ephot=0.5*alam**5*fphot
C         write(6,*) ' ephot = ', ephot, ' ebrems = ', ebrems
C         write(6,*) ' epair = ', epair
          fplas1 = (2.320d-7 + 8.449d-8*axi + 1.787d-8*axi**2)
     &     *exp(-.5646*axi)
          fplas2=axi**3 + 2.581d-2/alam + 1.734d-2/alam/alam
     &     + 6.990d-4/alam**3
          fplas=fplas1/fplas2
          eplas=rho*rho/8.d0*fplas 
C         write(6,*) ' eplas = ', eplas , ' eps = ', epsneu   
C     Munakata, Kohyama, and Itoh 1985 corrections
   
          eplas = eplas*CV2
          qphot=1.d0+rho/2.d0*(1.875d8*alam + 1.653d8*alam2
     & + 8.499d8*alam2*alam -1.604d8*alam4)**(-1.d0)
          qphot=0.666d0/qphot*(1.d0 + 2.045d0*alam)**(-2.066d0)
          fphot=fphot*0.893d9*alam4*alam4*(1.d0+143.8d0*alam**(3.555d0)
     &    )**(-0.3516d0)*axi**3*exp(0.556d0*axi**4.48d0/
     &    (150.d0+axi**3.3d0))
          fphot=fphot*0.5d0*((CV2+CA2))*(1.d0 - (CV2-CA2)/(CV2+CA2)
     &    *qphot)
          ephot=fphot/rho

C
          qpair=(10.7480d0*alam2 + 0.3967d0*sqrt(alam)
     &    +1.005d0)**(-1.d0)*(1.d0+rho/2.d0*(7.692d7*alam*alam2
     &    +9.715d6*sqrt(alam))**(-1.d0))**(-0.3d0)
          epair=epair*0.5d0*(CV2+CA2)*(1.d0+(CV2-CA2)/
     &    (CV2+CA2)*qpair) 
          
          

          epsneu=(eplas+ephot+ebrems+epair)

          E = E - epsneu
       endif

      if(xh.gt.1.d-37) then
c......1. Calculate screening factors fpp (weak&intermediate screening)
c.........a. Weak pp screening:

             kb=.5d0
             etahb=1.127d0
             b=1.d0
             zetab=2.d0
             mu=etahb**2/zeta
             Lambda0=1.88d+8*sqrt(Rho/(mu*T**3))
             H12=kb*etahb*zetab*(Lambda0**b)
             fppw=exp(H12)

c.........b. Weak he3he3 screening

             etahb=2.d0
             zetab=8.d0
             mu = etahb**2/zeta
             Lambda0=1.88d8*sqrt(Rho/(mu*T**3))
             h33= kb*etahb*zetab*(Lambda0**b)
             fhe3he3=exp(h33)
c             fhe3he3 = 1.1
c            if(H12.gt.(.1)) then
c..........c Intermediate pp screening:

                 b=0.860d0
                 kb=0.380d0

                 sumni=Z0
                 do i=1,14
                     n(i)=(xba(i)*Rho)/(HMASS)
                     sumni=sumni+n(i)
                 end do
                 avz3b=Z0
                 zbar=Z0
                 do i=1,14
                     avz3b=avz3b+Xz(i)**(3.d0*b-Z1)*n(i)
                     zbar=zbar+Xz(i)*n(i)
                 end do
                 avz3b=avz3b/sumni
                 zbar=zbar/sumni
                 etahb=avz3b/(1.127d0**(3.d0*b-2.d0)*zbar**(2.-2.*b))

                 zetab=1.63d0
                 H12=kb*etahb*zetab*(Lambda0**b)
                 fppi=exp(H12)
C
c            end if
          fpp = min(fppi,fppw)
C
c......2. Calculate effective cross section factors Spp and Spd, She3he3

       Ap=Z1
       Ad=2.d0
       Ahe3=3.d0
C
       N0=6.0225d+23
C
C       Np=Rho*(xh/Ap)*N0
C       nd=Rho*(xd/Ad)*N0

       App=.5d0
       Apd=(Ap*Ad)/(Ap+Ad)
       Ahe3he3 = Ahe3/2.d0
       S0pp=4.07d-22
       S0pd=2.5d-04
       S0he3he3 = 5.15d3

       dS0pp=4.52d-24
       dS0pd=7.9d-06
       dS0he3he3 = -0.9d0

       T6=T/1.d6
       taupp=42.487d0*(App/T6)**Z13
       taupd=42.487d0*(Apd/T6)**Z13
       tauhe3he3 =42.487d0*(16.d0*Ahe3he3/T6)**Z13

       E0pp=taupp*k*T*Z13
       E0pd=1.2204d0*(Apd*T6**2)**Z13
       E0he3he3 = tauhe3he3*k*T*Z13

       Spp=S0pp*(Z1+(5.d0/(12.d0*taupp)))+dS0pp*(E0pp+.97222d0*k*T)
       Spd=S0pd*(Z1+(5.d0/(12.d0*taupd)))+dS0pd*(E0pd+.97222d0*k*T)
       She3he3 = S0he3he3*(Z1+ (5.d0/(12.d0*tauhe3he3)))
     +   +dS0he3he3*(E0he3he3 + .97222d0*k*T)

c......Get the average products of cross section times velocity:

       sigvpp=1.3005d-15*(Z1/(App*T6**2))**Z13
       sigvpp=sigvpp*fpp*Spp*exp(-taupp)

       sigvpd=1.3005d-15*(Z1/(Apd*T6**2))**Z13
c      Xpd = Z0
       sigvpd=sigvpd*fpd*Spd*exp(-taupd)

       sigvhe3he3=1.3005d-15*(Zhe3**2/(Ahe3he3*T6**2))**Z13
       sigvhe3he3=sigvhe3he3*fhe3he3*She3he3*exp(-tauhe3he3)
c......Get the nuclear generation rates for the three reactions

       Rpp = np*np*sigvpp/2.d0
       Rhe3he3 = nhe3*nhe3*sigvhe3he3/2.d0
c......Correct nd for equilibrium deuterium burning
       Xpd=np*sigvpd
       if(Xpd*dtime.gt.100.d0) then
         nde=Z0
       else
         nde=nd*exp(-Xpd*dtime)
       endif
       if(Xpd*dtime.gt.0.01d0) then
         Rpd = (nd-nde)/dtime
       else
         Rpd = nde*Xpd
       endif

c......Get the Q values (adjusted for neutrino losses)

       Qpp = 2.310d-6 - 4.245d-7
       Qpd = 8.801d-6
       Qhe3he3 = 2.06d-5
c......Add on extra deuterium burning Q to the Qpp

       Qpp = Qpp + Qpd
c.......Add on he3 contribution to Qpp to form Qpp1
       Qpp1=Qhe3he3/2.d0 + Qpp

c......Get the energy generation rate

       if(iflip.eq.0) then
         E = E + Qpp*Rpp/Rho
         E = E + Qhe3he3*Rhe3he3/Rho
c......Add on contribution from primordial deuterium (Note that is is
c......implicitly assumed that deuterium is used up by the time Rpp is
c......significant.)
         E = E + Qpd*Rpd/Rho
         if(iwrite.gt.0)                                                &
     &     write(6,779) 'j2,nde,nd,Xpd,Rpd,E:',iwrite,nde,nd,Xpd,Rpd,E
 779       format(a,i5,1p,7E10.2)
       else

c....    He3 in statistical equilibrium.

c...     PPII and PPIII corrections due to Parker, Bahcall and Fowler, 1964.
C...     (Fit to curves for Y=0.1, 0.5, 0.9)
         xpp23=Z1+max(Z0,(T6-10.d0)/36.d0)
         xpp23=min(1.5d0,xpp23) + max(Z0,XHe4*(0.8d0-((T6-18.)/11.)**2))
         xpp23=min(1.95d0,xpp23)
c
         E = E + xpp23*Qpp1*Rpp/rho
c...     Compute CNO generation rate (from Bahcall 1989)
         Apc = 12./(12.+Z1)
         Apn = 14./(14.+Z1)
         Apo = 16./(16.+Z1)
         sigvpc=1.3005d-15*(6.d0/(Apc*T6**2))**Z13
         sigvpn=1.3005d-15*(7.d0/(Apn*T6**2))**Z13
         sigvpo=1.3005d-15*(8.d0/(Apo*T6**2))**Z13
         fpc=Z1
         fpn=Z1
         fpo=Z1
         taupc=136.93d0/T6**Z13
         taupn=152.31d0/T6**Z13
         taupo=166.96d0/T6**Z13
         S0pc=1.45
         S0pn=3.32
         S0po=9.4
         Spc=S0pc*(Z1+(5.d0/(12.d0*taupc)))
         Spn=S0pn*(Z1+(5.d0/(12.d0*taupn)))
         Spo=S0po*(Z1+(5.d0/(12.d0*taupo)))
         sigvpc=sigvpc*fpc*Spc*exp(-taupc)
         sigvpn=sigvpn*fpn*Spn*exp(-taupn)
         sigvpo=sigvpo*fpo*Spo*exp(-taupo)
         nc=xc*RHO/(HMASS*12.)
         nn=xn*RHO/(HMASS*14.)
         no=xo*RHO/(HMASS*16.)
         Qpc=1.7635d-5
         Qpn=2.2466d-5
         Qpo=5.6961d-6
         Xpc=np*sigvpc
         Xpn=np*sigvpn
         Xpo=np*sigvpo
         Rpc=nc*Xpc
         Rpn=nn*Xpn
         Rpo=no*Xpo
         Rcno=Rpc+Rpn+Rpo
         Ecno=(Qpc*Rpc+Qpn*Rpn+Qpo*Rpo)/RHO

         E = E + Ecno

       end if
      endif
C       write(6,1414)E 
C       write(6,1414)T,Rhe3he3,fppw,fppi,nP,rho,rpP
1414   format(1x,1p,8E11.3)
       if(ifirst.eq.1 .and. T.gt.2.e7) then
         ifirst=2
         write(6,'(" xc,xn,xo=",1p,3E12.4," Ecno=",E12.4," E=",E12.4)')
     &     xc,xn,xo,Ecno,E
        write(6,1414) T,Rhe3he3,fppw,fppi,nP,rho,Rpp
       endif
       return

       end

      subroutine compchange
      include 'parm.h'
      include 'var.h'
      dimension pp(5),AA(3,4)
c     dimension AA(3,4)
      COMMON/ABUND/XCA(14),H1(14),A(14)
      DATA Z0,Z1,Z13,HMASS/0.d0,1.d0,0.333333333333333333d0,1.6732D-24/
      data ifirst/0/
      save
c
C  Reset all abundances to original values when the deuterium
C  abundance is given as a negative value.
      if(XH2.lt.Z0) then
        H1(2)=abs(XH2)
        SUM = 1.
        DO J = 1,14
          SUM = SUM - H1(J)
        ENDDO
C H1(2) is temporarily used as the primordial deuterium abundance and
C then replaced by the initial helium abundance.
        H1(2)=SUM
        do i=1,N
          helium3(i) = 1.d-5
          hydrogen(i) = H1(1)+Z13*helium3(i)
          helium4(i) = H1(2)-4.d0/3.d0*helium3(i)
          deuterium(i)=abs(XH2)
          carbon(i)=H1(3)
          trogen(i)=H1(12)
          oxygen(i)=H1(13)
        enddo
        return
      endif
c
C  The following redistributes any inaccuracy of abundances to the
C  most abundant species.
c
      if(ifirst.eq.0) then
        ifirst=1
        XYZx = hydrogen(N) + helium4(N) + helium3(N) + deuterium(N)     &
     &       + oxygen(N) + trogen(N) + carbon(N)
        do i=1,N
          pp(1)=hydrogen(i)
          pp(2)=helium4(i)
          pp(3)=carbon(i)
          pp(4)=trogen(i)
          pp(5)=oxygen(i)
          XYZ= hydrogen(i) + helium4(i) + helium3(i) + deuterium(i)     &
     &       + oxygen(i) + trogen(i) + carbon(i)
          k=1
          do j=2,5
            if(pp(j).gt.pp(k)) k=j
          enddo
          pp(k)=pp(k) + XYZx-XYZ
          hydrogen(i) = pp(1)
          helium4(i)  = pp(2)
          carbon(i)   = pp(3)
          trogen(i)   = pp(4)
          oxygen(i)   = pp(5)
        enddo
      endif
      YYmax = hydrogen(N) + helium4(N) + helium3(N) + deuterium(N)
      do i=1,N 
        Tnuc = X(i,4)     
        olddeuterium=deuterium(i)
        oldhydrogen= hydrogen(i)
        oldhelium3 = helium3(i)
        oldhelium4 = helium4(i)
        oldcarbon  = carbon(i)
        oldtrogen  = trogen(i)
        oldoxygen  = oxygen(i)
        call flip(Tnuc,hydrogen(i),helium4(i),helium3(i),iflip)
        call invstate(41,x(i,1),x(i,4),hydrogen(i),helium4(i),Crad,RHO, &
     &               cP,alpha,beta,delta)
        QNh1 = hydrogen(i)*rho/HMASS
        QNh2 =deuterium(i)*rho/2.d0/HMASS
        QNhe3 = helium3(i)*rho/3.d0/HMASS

        call nucrat(Rho,Tnuc,QNh1,QNh2,QNhe3,dtime,ee,                  &
     &        Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                   &
     &        helium4(i),carbon(i),trogen(i),oxygen(i),                 &
     &        xpp23,iflip,0)

        if(Tnuc.gt.5.d7) then
          helium4(i)=helium4(i)+deuterium(i)+helium3(i)
          deuterium(i)=Z0
          helium3(i)=Z0
        else
          if(Xpd*dtime.gt.100.d0 .or. deuterium(i).lt.1.d-40) then
            deuterium(i)=Z0
          else
            deuterium(i)=deuterium(i)*exp(-Xpd*dtime)
          endif
          dxddt = (deuterium(i) - olddeuterium)/dtime
        endif
        dxhe4dt= -(3.d0*R3a+R12a)*4.d0*HMASS/RHO
        dxcdt  = (R3a-R12a)*12.d0*HMASS/RHO
        dxodt  = R12a*16.d0*HMASS/RHO
        if(iflip.eq.0) then
          dxhe3dt =(Rpp-2.d0*Rhe3)*3.d0*HMASS/rho-1.5d0*dxddt
          helium3(i) = dtime*dxhe3dt + helium3(i)
          if(helium3(i).lt.Z0) then
            helium3(i)=Z0
            dxhdt = (-3.d0*Rpp)*HMASS/rho
            dxhe4dt= dxhe4dt - dxhdt 
            helium4(i) = helium4(i) + oldhelium3
            dxh=min(oldhelium*Z13,hydrogen(i))
            hydrogen(i)= hydrogen(i) - dxh
            helium4 (i)= helium4 (i) + dxh
          else
            dxhdt  = (-3.d0*Rpp+2.d0*Rhe3)*HMASS/rho+.5d0*dxddt
            dxhe4dt= dxhe4dt + Rhe3 *4.d0*HMASS/rho
          endif
        else
          brem=Z1/dtime
          AA(1,1)=brem+Xpc
          AA(1,2)=-(Z1-4.d-4)*Xpn
          AA(1,3)=Z0
          AA(1,4)=vcarbon(i)*brem
          AA(2,1)=-Xpc
          AA(2,2)=brem+Xpn
          AA(2,3)=-Xpo
          AA(2,4)=vtrogen(i)*brem
          AA(3,1)=Z0
          AA(3,2)=-4.d-4*Xpn
          AA(3,3)=brem+Xpo
          AA(3,4)=voxygen(i)*brem
          CALL GIRL(AA,3,1)
          carbon(i)=AA(1,4)
          trogen(i)=AA(2,4)
          oxygen(i)=AA(3,4)
          dxhdt  = -(2.d0*xpp23*Rpp+4.d0*Rcno)*HMASS/rho
c         if(i.eq.1) then
c           write(6,*) xpp23,Rpp,HMASS,rho,Rcno,-dxhdt
c           stop
c         endif
          dxhe4dt= dxhe4dt - dxhdt 
        end if

        hydrogen(i)= dtime*dxhdt   + hydrogen(i)
        helium4(i) = dtime*dxhe4dt + helium4(i)
        if(oldhydrogen.gt.Z0 .and. hydrogen(i) .lt. 1.d-37) then
          hydrogen(i) = Z0
          helium3(i) = Z0
          helium4(i) = YYmax
        endif
        if(oldhelium4.gt.Z0 .and. helium4(i).lt.1.d-37) then
          helium4(i)= Z0
          carbon(i) = carbon(i) + oldhelium4
        else
          carbon(i) = carbon(i) + dxcdt*Dtime
        endif
        if(oldcarbon.gt.Z0 .and. carbon(i).lt.Z0) then
          carbon(i)= Z0
          oxygen(i) = oxygen(i) + carbon(i)-dxcdt*Dtime
        else
          oxygen(i) = oxygen(i) + dxodt*Dtime
        endif
c       if(i.eq.1) then
c         Znet=1.d0-oldhydrogen-oldhelium4-oldhelium3-olddeuterium      &
c    &      -oldcarbon-oldnigrogen-oldoxygen
c         write (6,200) i,oldhydrogen,oldhelium4,oldhelium3,            &
c    &        olddeuterium,oldcarbon,oldnigrogen,oldoxygen,Znet,iflip
c         Znet=1.d0-hydrogen(1)-helium4(1)-helium3(1)-deuterium(1)      &
c    &      -carbon(1)-trogen(1)-oxygen(1)
c         write (6,200) i,hydrogen(1),helium4(1),helium3(1),            &
c    &        deuterium(1),carbon(1),trogen(1),oxygen(1),Znet
c200      format(i4,' H,He4,He3,D,C,N:=',8F8.5,i2)
c       endif
      end do 
 
C..... MixinG

      iset=0
      ic1=0
      ic2=0
      do j=1,N
        if(iconvect(j) .eq. 0) then
          if(iset.eq.0) then
            ic1=j
            iset=1
          else
            ic2=j
          end if
        else
          if(iset.eq.1) then
            ic2=j-1
            iset=2
          end if
        end if
        if(iset.eq.2 .or. (iset.eq.1 .and. j.eq.N)) then
          if(ic1.lt.ic2) then
            runsum = Z0
            runsum1= Z0
            runsum2= Z0
            runsum3= Z0
            runsum4= Z0
            runsum5= Z0
            runsum6= Z0
            runsum7= Z0
            if(ic2 .eq. N .and. Zflux.gt.0.) then
              runsum  = Zflux*dtime
              runsum1 = H1(1)*runsum
              runsum2 = XH2 *runsum
              runsum3 = 1.d-5*runsum
              runsum4 = H1(2)*runsum
              runsum5 = H1(3)*runsum
              runsum6 = H1(12)*runsum
              runsum7 = H1(13)*runsum
            endif
            do i = ic1,ic2
              runsum1 = runsum1 + hydrogen(i)*dM(I) 
              runsum2 = runsum2 + deuterium(i)*dM(I) 
              runsum3 = runsum3 + helium3(i)*dM(I) 
              runsum4 = runsum4 + helium4(i)*dM(I)
              runsum5 = runsum5 + carbon(i)*dM(I)
              runsum6 = runsum6 + trogen(i)*dM(I)
              runsum7 = runsum7 + oxygen(i)*dM(I)
              runsum  = runsum  + dM(i)
            enddo
            Xh = runsum1/runsum
            xd = runsum2/runsum
            xhe3=runsum3/runsum
            xhe4=runsum4/runsum
            xc  =runsum5/runsum
            xn  =runsum6/runsum
            xo  =runsum7/runsum
C            write(6,*) xh,xd,xhe3,xhe4
            do i = ic1,ic2
              hydrogen(i) = xh
              deuterium(i)= xd
              helium3(i)  = xhe3
              helium4(i)  = xhe4
              carbon(i)   = xc
              trogen(i)   = xn
              oxygen(i)   = xo
            end do
          endif
          iset=0
          ic1=0
          ic2=0
        endif
      enddo
c     Znet=1.d0-hydrogen(1)-helium4(1)-helium3(1)-deuterium(1)          &
c    &      -carbon(1)-trogen(1)-oxygen(1)
c     write (6,200) hydrogen(1),helium4(1),helium3(1),deuterium(1),     &
c    &        carbon(1),trogen(1),oxygen(1),Znet
c
C     END MIXING 
      RETURN
      END 
      subroutine flip(T,hydrogen,helium4,helium3,iflip)
      include 'parm.h'
      data Z0/0.d0/
      save
c
      if(T.gt.1.0d+7 .or. helium3.eq.Z0) then
        iflip=1
        helium4=helium4+helium3
        dX=min(helium3/3.d0,hydrogen)
        hydrogen=hydrogen-dX
        helium4 =helium4 +dX
        helium3=Z0
      else
        iflip=0
      end if
      return
      end
