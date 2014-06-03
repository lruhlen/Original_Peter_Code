      program polytr
      IMPLICIT real*8 (A-H,O-Z), INTEGER*4(I-N)
      PARAMETER (MJ=3000,MH=4)
      dimension zM(MJ),R(MJ),zL(MJ),T(MJ),RHO(MJ)
      data PI4,GRAV/12.566371,6.6704E-8/
      data Rsun,Rgas/6.96E10,82985077./
      data DmassF/1.d-5/
c
      open(3,file='polyout',status='unknown')
      N=MJ
      call comrd
      read (5,301) Zmas,Rstar,Teff,XX,YY,Pexp
c     write(6,201) Zmas,Rstar,Teff,XX,YY,Pexp
 301  format(7x,1p,e9.2,7x,e9.2,6x,e9.2,3x,0p,f9.6,3x,f9.6,3x,f6.3)
 201  format(' Zmass=',1p,e9.2,' Rstar=',e9.2,' Teff=',e9.2,            &
     & ' X=',0p,f9.6,' Y=',f9.6,' E=',f7.4)
      amu = 2. * XX + 0.75 * YY + 0.5 * (1.-XX-YY)
      Rmu=Rgas*amu
      Zfac=1.0
      Zmass=Zmas
      kmax=1
      do kkk=1,kmax
        Zmass=Zmas*Zfac
        A=3.65375/Rstar
        RHOc=5.99071*3.*Zmass/(PI4*Rstar**3)
        zK=PI4*GRAV/(1.+Pexp)*RHOc**(1.-1./Pexp)/A**2
        Pc=zK*RHOc**(1.+1./Pexp)
        write(6,202) RHOc,Pc,A,zK
 202    format(' RHOc=',1p,e11.4,' Pc=',e11.4,' A=',e11.4,' K=',e11.4)
        Dmass = zMass*DmassF
c       RR0   = 6.*Pc*(1.+1./Pexp)/(PI4*GRAV*RHOc**2)
c       RR1   = ((1.5*Dmass)/(PI4*RHOc**2))**(2./3.)
c       RHO(1)= RHOc*(1.-RR1/RR0)
c       P     = Pc*(1.-(1.+1./Pexp)*RR1/RR0)
        RHO(1)= RHOc
        P     = Pc
        T(1)  = P/(RHO(1)*Rmu)
        zL(1) = Dmass * T(1) 
        R3    = (3.*Dmass/(PI4*RHO(1)))
        zM(1) = Dmass
        R(1)  = R3**(1./3.)
        DD    = 0.
        do j=2,N
          Dmold=Dmass
          if(DD.gt.0.05) Dmass=max(Dmass*0.9,DmassF*Zmass)
          if(DD.lt.0.04) Dmass=min(Dmass*1.05,0.01*Zmass)
          if(DD.lt.0.02) Dmass=min(Dmass*1.5,0.01*Zmass)
          Dmhf  = 0.5*(Dmold+Dmass)
          dP    = -GRAV*zM(j-1)*Dmhf/(PI4*R(j-1)**4)
          if(P+dP .le. 0.) goto 998
          P     = P+dP
          zM(j) = zM(j-1)+Dmass
          RHO(j)= (P/zK)**(Pexp/(1.+Pexp))
          T(j)= P/(RHO(j)*Rmu)
          zL(j) =zL(j-1) + Dmass*T(j)
          dR3=3.*Dmass/(PI4*RHO(j))
          R3    = R3+dR3
          R(j) = (R(j-1)**3 + dR3)**(1./3.)
          if(kkk.eq.kmax) write(6,203) J,zM(j),R(j),RHO(j),P
 203      format(1p,i5,5E12.4)
          if(T(j).lt.7.0E4) go to 999
          DD    = abs(dP/(P+dP)) + abs(dR3/R3)/3.
        enddo
        goto 999
 998    j=j-1
 999    continue
        XL = PI4*5.67e-5*Teff**4*Rstar**2 
        Zfac=Zmas/zM(J)
        write(6,203) J,Zmass,Zfac,XL
      enddo
      write(6,*) Teff,XL,zL(J),N
      FAC=zL(J)/XL 
      do i = 1,J
        zL(i)= zL(i)/fac/3.85e33
      end do 
      write(3,204) (zM(i)/zM(J),R(i)/Rsun,T(i),zL(i),RHO(i),i=1,J)
 204  format(1p,5E12.4) 
      write(3,204) Zmas
      stop 
      end
      subroutine comrd
      character*1 line(1)
      character*80 lin
      equivalence (lin,line)
c
1     read (5,203,end=999) lin
      if(line(1).eq.' ') then
      write(6,203)         lin
        backspace 5
        return
      else
        if(line(1).ne.'c' .and. line(1).ne.'C')
     *  write(6,203)  lin
c202    format(1x,a)
      endif
      goto 1
 999  write(6,203) 'COMRD: END OF INPUT DATA'
 203  format(a)
      return
      end
