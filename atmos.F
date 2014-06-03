      SUBROUTINE ATMOS(Tatm,RHOatm,Ratm,Patm,Rstar,Xlum,jwrite)
      include 'parm.h'
      include 'var.h'
      parameter  (MTAU=3000)
      CHARACTER*1 CS(2)
      DIMENSION TAU(MTAU),TTAU(MTAU),RHOTAU(MTAU),PTAU(MTAU),           &
     & RTAU(MTAU),ZMTAU(MTAU),ZATG(MTAU),ZRADG(MTAU),                   &
     & ZTRUG(MTAU),ICV(MTAU)
c
c     DATA Arad3/2.52197145D-15/
      DATA SIG,RG,PI/5.67051d-5,8.31451d7,3.14159265359d0/
      DATA GRAV,CC/6.6704D-8,2.99792458d10/
      DATA Z0,Z05,Z1,Z2,Z10,Z13,Z23/0.d0,0.5d0,1.d0,2.d0,10.d0,         &
     & .333333333333333333d0,.6666666666666666667d0/
      DATA CS/'*',' '/
      data ifirst/0/
      SAVE
c
      iwrite=jwrite
      kwrite=abs(iwrite)
      ATMASS1=0.90d0*dM(N)
      ATMASS2=dM(N)
      Arad3=SIG*4.d0/(CC*3.d0)
      XX=hydrogen(N)
      YY=helium3(N)+helium4(N)
      DELTAU = .001d0
      dTAU05 = Z05*DELTAU
      ITAU23 = 0
      IDELM  = 0
      RTAU23 = Rstar
      ZTAU23 = zM(N)
      TEFF4  = Xlum/(4.d0*PI*SIG*Rstar**2)
      TEFF   = SQRT(SQRT(TEFF4))
      RK0    = Z0
      ZK0    = Z0
      TK0    = SQRT (SQRT(TEFF4*(0.75d0*dTAU05 + Z05)))
      Rat    = Rstar
      Zat    = zM(N)
C
C   Define atmosphere values for JK=1 at 1/2 DELTAU...
C
      Prad=Crad*Arad3*TK0**4
      Tlog=log10(TK0)
      Pout=max(Prad*Z2,Prad+10.)
      call invstate(11,Pout,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
      RHOlog = log10(RHO)
      call opacity(N,Tlog,RHOlog,bkap) 
      AKK = Z10**bkap
      GD  = DELTAU*GRAV*Zat/Rat**2
      PK0 = Pout + GD/AKK
      call invstate(12,PK0,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
      RHOlog = log10(RHO)
      call opacity(N,Tlog,RHOlog,bkap) 
      AKK = Z10**bkap
      do itr=1,30
        G0 = (PK0-Pout)*AKK - GD
        PK1=PK0*1.001d0
        call invstate(13,PK1,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
        G1 = (PK1-Pout)*AKK - GD
        dGdP=(G0-G1)/(PK0-PK1)
        DELP = -G0/dGdP
        DELP = max(DELP,-Z05*PK0)
        DELP = min(DELP,.9d0*PK0)
        DELP = max(DELP,.8d0*(Prad-PK0))
        if(itr.gt.11) DELP=Z05*DELP
        PK0 = PK0 + DELP
        call invstate(14,PK0,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
        if(itr.gt.2 .and. abs(DELP)/(PK0+PK1) .lt. 1.e-5) goto 3
      enddo
 3    continue
      if(ifirst.eq.0 .or. itr.ge.30)
     &  write(6,'(a,1p,4e12.4,a,i3,a)')                                 &
     &     'ATMOS: First values of P,RHO,T,AKM:',PK0,RHO,TK0,AKK,       &
     &     ' after ',itr,' iterations'
      ifirst=1
      TAU(1)   = dTAU05
      PTAU(1)  = PK0
      RTAU(1)  = Z0
      ZMTAU(1) = Z0
      TTAU(1)  = TK0
      RHOTAU(1)= RHO
      call nabla(PK0,Rat,Xlum,TK0,Zat,TPNAB,ADNAB,RADNAB,N,IC)
      ZATG(1)  = ADNAB
      ZRADG(1) = RADNAB
      ZTRUG(1) = TPNAB*PK0/TK0
      ICV(1)   = IC
C
C   Done defining atmospheric values at outermost point. Now do the
C   rest of the atmosphere using Runge-Kutta integration...
C
      do JK=2,MTAU
C
C   First Runge-Kutta step: get K1 values for P (pressure), R (radius),
C   Z (mass), and T (temperature) by straightfoward extrapolation with
C   a full step of DELTAU...  R and Z are calculated from TAU = 0
C
        TAU(JK) = TAU(JK-1) + DELTAU
C
C   The following two lines determine radius and mass from center of
C   star, except when outside TAU = 2/3, where Rstar and Mstar are 
C   used.
C
        Rat = Rstar - max(Z0,RK0-RTAU23)
        Zat = zM(N) - max(Z0,ZK0-ZTAU23)
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK1 = GD
        dRK1 = dTAK/RHO
        dZK1 = dTAK*4.d0*PI*Rat**2
        dTK1 = GD*TPNAB
c       if(JK.eq.2)
c    &  dTK1 = SQRT (SQRT(TEFF4*(0.75d0*TAU(JK) + Z05))) - TK0
C
C   Second Runge-Kutta step: get K2 values for P, R, Z and T at 1/2 step
C
        PM2 = PK0 + Z05*dPK1
        RM2 = RK0 + Z05*dRK1
        ZM2 = ZK0 + Z05*dZK1
        TM2 = TK0 + Z05*dTK1
        Rat = Rstar - max(Z0,RM2-RTAU23)
        Zat = zM(N) - max(Z0,ZM2-ZTAU23)
        call nabla(PM2,Rat,Xlum,TM2,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        call invstate(15,PM2,TM2,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TM2)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK2 = GD
        dRK2 = dTAK/RHO
        dZK2 = dTAK*4.d0*PI*Rat**2
        dTK2 = GD*TPNAB
C
C   Third Runge-Kutta step: get K3 values for P, R, Z and T at 1/2 step
C
        PM2 = PK0 + Z05*dPK2
        RM2 = RK0 + Z05*dRK2
        ZM2 = ZK0 + Z05*dZK2
        TM2 = TK0 + Z05*dTK2
        Rat = Rstar - max(Z0,RM2-RTAU23)
        Zat = zM(N) - max(Z0,ZM2-ZTAU23)
        call nabla(PM2,Rat,Xlum,TM2,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        call invstate(16,PM2,TM2,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TM2)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK  = Z10**bkap
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK3 = GD
        dRK3 = dTAK/RHO
        dZK3 = dTAK*4.d0*PI*Rat**2
        dTK3 = GD*TPNAB
C
C   Fourth Runge-Kutta step: get K4 values for P, R, Z and T at full step
C
        PM2 = PK0 + dPK3
        RM2 = RK0 + dRK3
        ZM2 = ZK0 + dZK3
        TM2 = TK0 + dTK3
        Rat = Rstar - max(Z0,RM2-RTAU23)
        Zat = zM(N) - max(Z0,ZM2-ZTAU23)
        call nabla(PM2,Rat,Xlum,TM2,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        call invstate(17,PM2,TM2,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TM2)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK  = Z10**bkap
        dTAK = DELTAU/AKK
        GD   = dTAK*GRAV*Zat/Rat**2
        dPK4 = GD
        dRK4 = dTAK/RHO
        dZK4 = dTAK*4.d0*PI*Rat**2
        dTK4 = GD*TPNAB
C
C   Final Runge-Kutta step: Combine K1, K2, K3, and K4 values for P, R, Z
C   and T at the grid point JK
C
        PK0 = PK0 + (dPK1 + Z2*(dPK2 + dPK3) + dPK4)/6.d0
        RK0 = RK0 + (dRK1 + Z2*(dRK2 + dRK3) + dRK4)/6.d0
        ZK0 = ZK0 + (dZK1 + Z2*(dZK2 + dZK3) + dZK4)/6.d0
        TK0 = TK0 + (dTK1 + Z2*(dTK2 + dTK3) + dTK4)/6.d0
        Rat = Rstar - max(Z0,RK0-RTAU23)
        Zat = zM(N) - max(Z0,ZK0-ZTAU23)
        call invstate(18,PK0,TK0,XX,YY,Crad,RHO,cP,alpha,Rbeta1,delta)
        Tlog   = log10(TK0)
        RHOlog = log10(RHO)
        call opacity(N,Tlog,RHOlog,bkap) 
        AKK = Z10**bkap
        call nabla(PK0,Rat,Xlum,TK0,Zat,TPNAB,ADNAB,RADNAB,N,IC)
        PTAU(JK)  = PK0
        RTAU(JK)  = RK0
        ZMTAU(JK) = ZK0
        TTAU(JK)  = TK0
        RHOTAU(JK)= RHO
        ZATG(JK)  = ADNAB
        ZRADG(JK) = RADNAB
        ZTRUG(JK) = TPNAB*PK0/TK0
        ICV(JK)   = IC
        if(TAU(JK).ge.Z23 .and. ITAU23.eq.0) then
          ITAU23=JK
          FAC    = (TAU(JK) - Z23)/(TAU(JK)-TAU(JK-1))
          RTAU23 = FAC*RTAU(JK-1) + (Z1-FAC)*RTAU(JK)
          ZTAU23 = FAC*ZMTAU(JK-1) + (Z1-FAC)*ZMTAU(JK)
        endif
        if(ZK0 .GE. ATMASS1 .and. IDELM.eq.0) THEN
          IDELM=JK
          FAC = (ZK0 - ATMASS1)/(ZK0-ZMTAU(JK-1))
          Tatm   = FAC*TTAU  (JK-1) + (Z1-FAC)*TTAU  (JK)
          RHOatm = FAC*RHOTAU(JK-1) + (Z1-FAC)*RHOTAU(JK)
          Patm   = FAC*PTAU  (JK-1) + (Z1-FAC)*PTAU  (JK)
        endif
        if(ZK0 .GE. ATMASS2) GO TO 75
        IF(TAU(JK) .GT. Z10) THEN
          DELTA = min(1.08d0 * DELTAU , dTAUmx)
          DELTAU = max(DELTA , 1.01d0* DELTAU )
          DELTAU = 1.01d0* DELTAU
        ELSE
          IF(TAU(JK) .GT. .1d0) DELTAU = .01d0
          IF(TAU(JK) .GT. .8d0) DELTAU = .05d0
        ENDIF
      enddo
      goto 998
75    CONTINUE
      FAC = (ZK0 - ATMASS2)/(ZK0-ZMTAU(JK-1))
      Ratm = FAC*RTAU(JK-1)+ (Z1-FAC)*RTAU(JK) - RTAU23
      if(iwrite.ne.0) then
        write(6,300) Tatm,Ratm,Patm,RHOatm,TEFF,Rstar,Xlum
 300    format(/,' Atmospheric Parameters:  T=',1p,E10.2,               &
     &         ' R=',E10.2,' P=',E10.2,' RHO=',E10.2,/,                 &
     &         '     Stellar Parameters:  T=',E10.2,' R=',E10.2,        &
     &         ' L=',E10.2,/                                            &
     & '    J C    dM        M        P        R       dR        T',    &
     & '       RHO      TAU      ATG      RAD      TRU    1-beta')
        do j=JK,1,-1
          if(mod(j,kwrite).eq.0 .or. j.ge.JK-1 .or. j.le.2 .or.         &
     &      j.eq.IDELM  .or. j.eq.IDELM -1 .or.                         &
     &      j.eq.ITAU23 .or. j.eq.ITAU23-1) then
            Prad=Crad*Arad3*TTAU(j)**4
            beta=Prad/PTAU(j)
            write(6,301) j,CS(1+ICV(j)),ZMTAU(j),zM(N)-ZMTAU(j),PTAU(j),&
     &                 Rstar-RTAU(j)+RTAU23,RTAU(j),TTAU(j),RHOTAU(j),  &
     &                 TAU(j),ZATG(j),ZRADG(j),ZTRUG(j),beta
 301        format(1x,i4,1x,a1,1p,12E9.2)
 302        format(1x,4('*'),1x,a1,1p,12E9.2)
          endif
          if(j.eq.JK) then
            FAC = (ZMTAU(j) - ATMASS2)/(ZMTAU(j)-ZMTAU(j-1))
            PX  = FAC*PTAU(j-1) + (Z1-FAC)*PTAU(j)
            TX  = FAC*TTAU(j-1) + (Z1-FAC)*TTAU(j)
            RX  = FAC*TAU (j-1) + (Z1-FAC)*TAU (j)
            ZMX = FAC*ZMTAU(j-1) + (Z1-FAC)*ZMTAU(j)
            RHX = FAC*RHOTAU(j-1) + (Z1-FAC)*RHOTAU(j)
            ATGX= FAC*ZATG(j-1) + (Z1-FAC)*ZATG(j)
            RADX= FAC*ZRADG(j-1) + (Z1-FAC)*ZRADG(j)
            TRUX= FAC*ZTRUG(j-1) + (Z1-FAC)*ZTRUG(j)
            IC  = FAC*ICV(j-1) + (Z1-FAC)*ICV(j)
            beta= Crad*Arad3*(FAC/PTAU(j-1)*TTAU(j-1)**4                &
     &          +        (Z1-FAC)/PTAU(j  )*TTAU(j  )**4)
            write(6,302) CS(1+IC),ZMX,zM(N)-ZMX,PX,                     &
     &            Rstar-Ratm,Ratm+RTAU23,TX,RHX,RX,ATGX,RADX,TRUX,beta
          endif
          if(j.eq.IDELM) then
            FAC = (ZMTAU(j) - ATMASS1)/(ZMTAU(j)-ZMTAU(j-1))
            RX  = FAC*RTAU(j-1) + (Z1-FAC)*RTAU(j)
            TX  = FAC*TAU (j-1) + (Z1-FAC)*TAU (j)
            ZMX = FAC*ZMTAU(j-1) + (Z1-FAC)*ZMTAU(j)
            ATGX= FAC*ZATG(j-1) + (Z1-FAC)*ZATG(j)
            RADX= FAC*ZRADG(j-1) + (Z1-FAC)*ZRADG(j)
            TRUX= FAC*ZTRUG(j-1) + (Z1-FAC)*ZTRUG(j)
            IC  = FAC*ICV(j-1) + (Z1-FAC)*ICV(j)
            beta= Crad*Arad3*(FAC/PTAU(j-1)*TTAU(j-1)**4                &
     &          +        (Z1-FAC)/PTAU(j  )*TTAU(j  )**4)
            write(6,302) CS(1+IC),ZMX,ZM(N)-ZMX,Patm,                   &
     &            Rstar-RX+RTAU23,RX,Tatm,RHOatm,TX,ATGX,RADX,TRUX,beta
          endif
          if(j.eq.ITAU23) then
            FAC = (TAU(j) - Z23)/(TAU(j)-TAU(j-1))
            PX  = FAC*PTAU(j-1) + (Z1-FAC)*PTAU(j)
            TX  = FAC*TTAU(j-1) + (Z1-FAC)*TTAU(j)
            ZMX = FAC*ZMTAU(j-1) + (Z1-FAC)*ZMTAU(j)
            RHX = FAC*RHOTAU(j-1) + (Z1-FAC)*RHOTAU(j)
            ATGX= FAC*ZATG(j-1) + (Z1-FAC)*ZATG(j)
            RADX= FAC*ZRADG(j-1) + (Z1-FAC)*ZRADG(j)
            TRUX= FAC*ZTRUG(j-1) + (Z1-FAC)*ZTRUG(j)
            IC  = min(ICV(j-1),ICV(j))
            beta= Crad*Arad3*(FAC/PTAU(j-1)*TTAU(j-1)**4                &
     &          +        (Z1-FAC)/PTAU(j  )*TTAU(j  )**4)
            write(6,302) CS(1+IC),ZMX,zM(N)-ZMX,PX,                     &
     &                 Rstar,RTAU23,TX,RHX,Z23,ATGX,RADX,TRUX,beta
          endif
        enddo
      endif
      taujk=TAU(JK)
      RETURN
 998  write(6,*) ' atmos: 998'
      write(6,301) JK,CS(2),dTAUmx,TAU(JK-1),TAU(JK-2)
      write(6,*) Tatm,RHOatm,Ratm,Patm,Rstar,Xlum,taujk,iwrite
      stop 'ATMOS: 998'
      END
