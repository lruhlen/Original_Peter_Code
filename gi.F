      subroutine GI(K2)
c
      include 'parm.h'
      include 'var.h'
      DATA GRAV/6.6704D-8/
      DATA PI  /3.14159265359d0/ 
c     DATA ELOG/ .434294d0/
      DATA XINC/.0005d0/
      DATA HMASS/1.6732d-24/
      save
C
      J2=abs(K2)
      J1=J2-1
      J3=J2+1
c     if(vx(j2,4).gt.5.e6) then
c       I1=-J1
c       I2=-J2
c       I3=-J3
c     else
c       I1=J1
c       I2=J2
c       I3=J3
c     endif
      if(K2.gt.0) then
        do j=1,NG
        do k=1,NG
          HC(j,k)=0.d0
          HD(j,k)=0.d0
          HE(j,k)=0.d0
        enddo
        enddo
      endif
      if(dTIME.le.0.d0) then
        brem=0.d0
      else
        brem=Cwrk/dTIME
      endif
      if (J2.eq.1) then
        zMhf2=zM(J2)-0.5*dM(J2)
        zMhf3=zM(J2)+0.5*dM(J3)
        GMdM = GRAV*zM(J2)*(dM(J2)+dM(J3))/(8.d0*PI)
        call invstate(20,x(j3,1),x(j3,4),hydrogen(j3),helium4(j3),Crad, &
     &                RHO3,cP3,alpha3,beta3,delta3)
        call invstate(21,x(j2,1),x(j2,4),hydrogen(j2),helium4(j2),Crad, &
     &                RHO2,cP,alpha,beta,delta)
        call flip(x(j2,4),hydrogen(j2),helium4(j2),helium3(j2),iflip)
        QNh1 = hydrogen(j2)*RHO2/HMASS
        QNh2 =deuterium(j2)*RHO2/2.d0/HMASS
        QNhe3 = helium3(j2)*RHO2/3.d0/HMASS
        CALL NUCRAT(RHO2,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EP,              &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
        call nabla(x(j3,1),X(J3,2),X(J3,3),X(J3,4),zMhf3,               &
     &     TPNAB3,ADNAB,RADNAB,J3,IC)
        iconvect(j3) = ic 
        call nabla(X(J2,1),X(J2,2),X(J2,3),X(J2,4),zMhf2,               &
     &     TPNAB2,ADNAB,RADNAB,J2,IC)
        iconvect(j2) = ic 
        if(K2.gt.0) then
          PN=X(J2,1)*(1.d0+XINC)
          PM=X(J3,1)*(1.d0+XINC)
          call invstate(23,PN,x(j2,4),hydrogen(j2),helium4(j2),Crad,    &
     &                 RHP2,cPP,alphaP,beta,deltaP)
          call nabla(PN,X(J2,2),X(J2,3),X(J2,4),zMhf2,                  &
     &     TPNAB2P,ADNAB,RADNAB,J2,IC)
          call nabla(PM,X(J3,2),X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3P,ADNAB,RADNAB,J3,IC)
          RN=X(J2,2)*(1.d0+XINC)
          RM=X(J3,2)*(1.d0+XINC)
          call nabla(X(J2,1),RN,X(J2,3),X(J2,4),zMhf2,                  &
     &     TPNAB2R,ADNAB,RADNAB,J2,IC)
          call nabla(X(J3,1),RM,X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3R,ADNAB,RADNAB,J3,IC)
          XLN=X(J2,3)*(1.d0+XINC)
          XLM=X(J3,3)*(1.d0+XINC)
          call nabla(X(J2,1),X(J2,2),XLN,X(J2,4),zMhf2,                 &
     &     TPNAB2L,ADNAB,RADNAB,J2,IC)
          call nabla(X(J3,1),X(J3,2),XLM,X(J3,4),zMhf3,                 &
     &     TPNAB3L,ADNAB,RADNAB,J3,IC)
          TN=X(J2,4)*(1.d0-XINC)
          TM=X(J3,4)*(1.d0-XINC)
          call invstate(25,x(j2,1),TN,hydrogen(j2),helium4(j2),Crad,    &
     &                  RHT2,cPT,alphaT,beta,deltaT)
          CALL NUCRAT(RHO2,TN,QNh1,QNh2,QNhe3,DTIME,EPT,                &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          call nabla(X(J2,1),X(J2,2),X(J2,3),TN,zMhf2,                  &
     &     TPNAB2T,ADNAB,RADNAB,J2,IC)
          call nabla(X(J3,1),X(J3,2),X(J3,3),TM,zMhf3,                  &
     &     TPNAB3T,ADNAB,RADNAB,J3,IC)
          RHD=RHO2*(1.d0+XINC)
          QNh1 = hydrogen(j2)*rhD /HMASS
          QNh2 =deuterium(j2)*rhD /2.d0/HMASS
          QNhe3 = helium3(j2)*rhD /3.d0/HMASS
          CALL NUCRAT(RHD,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EPD,            &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          dRHdP2=(RHP2-RHO2)/(XINC*X(J2,1))
          dRHdT2=(RHT2-RHO2)/(-XINC*X(J2,4))
          dEPdR=(EPD-EP)/(XINC*RHO2)
          dEPdT=(EPT-EP)/(-XINC*X(J2,4)) + dEPdR*dRHdT2
          dEPdP=dEPdR*dRHdP2
          dcPdP=(cPP-cP)/(XINC*X(J2,1))
          dcPdT=(cPT-cP)/(-XINC*X(J2,4))
          ddeldP=(deltaP-delta)/(XINC*X(J2,1))
          ddeldT=(deltaT-delta)/(-XINC*X(J2,4))
          dNAB2P=(TPNAB2P-TPNAB2)/(XINC*X(J2,1))
          dNAB2R=(TPNAB2R-TPNAB2)/(XINC*X(J2,2))
          dNAB2L=(TPNAB2L-TPNAB2)/(XINC*X(J2,3))
          dNAB2T=(TPNAB2T-TPNAB2)/(-XINC*X(J2,4))
          dNAB3P=(TPNAB3P-TPNAB3)/(XINC*X(J3,1))
          dNAB3R=(TPNAB3R-TPNAB3)/(XINC*X(J3,2))
          dNAB3L=(TPNAB3L-TPNAB3)/(XINC*X(J3,3))
          dNAB3T=(TPNAB3T-TPNAB3)/(-XINC*X(J3,4))
          HD(1,1) =-1.d0
          HD(1,2) =-4.d0*GMdM/X(J2,2)**5
          HE(1,1) = 1.d0
          HD(2,1) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdP2
          HD(2,2) = 3.d0*X(J2,2)**2
          HD(2,4) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdT2
          HD(3,1) =-dM(J2)*( dEPdP                                      &
     &             -dcPdP * (X(J2,4)-VX(J2,4))*brem                     &
     &             +delta/RHO2 *brem                                    &
     &             +ddeldP/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdP2 * (X(J2,1)-VX(J2,1))*brem )
          HD(3,3) = 1.d0
          HD(3,4) =-dM(J2)*( dEPdT                                      &
     &             -cP * brem                                           &
     &             -dcPdT * (X(J2,4)-VX(J2,4))*brem                     &
     &             +ddeldT/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdT2 * (X(J2,1)-VX(J2,1))*brem )
          HD(4,1) = 0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB2P
          HD(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2R
          HD(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2L
          HD(4,4) =-1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB2T
          HE(4,1) =-0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB3P
          HE(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3R
          HE(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3L
          HE(4,4) = 1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB3T
        endif
        G(1) = X(J3,1)-X(J2,1) + GMdM/X(J2,2)**4
        G(2) = X(J2,2)**3 - dM(J2)*0.75d0/(PI*RHO2)
        G(3) = X(J2,3) - dM(J2)*( EP - cP * (X(J2,4)-VX(J2,4))*brem     &
     &         + delta/RHO2 * (X(J2,1)-VX(J2,1))*brem )
        G(4) = (X(J3,4)-X(J2,4))                                        &
     &         -(X(J3,1)-X(J2,1))*0.5d0*(TPNAB3+TPNAB2)
C       write(6,*) j2, ' Gs ', g(1),g(2),g(3),g(4) 
        return
      endif 
      if(J2.gt.1 .and. J2.lt.N) then
        zMhf2=zM(J2)-0.5*dM(J2)
        zMhf3=zM(J2)+0.5*dM(J3)
        GMdM = GRAV*zM(J2)*(dM(J2)+dM(J3))/(8.d0*PI)
        if(K2.gt.0) then
          TPNAB2= TPNAB3
        else
          call nabla(X(J2,1),X(J2,2),X(J2,3),X(J2,4),zMhf2,             &
     &       TPNAB2,ADNAB,RADNAB,J2,IC)
        endif
        call invstate(26,x(j2,1),x(j2,4),hydrogen(j2),helium4(j2),      &
     &                Crad,RHO2,cP,alpha,beta,delta)
        call flip(x(j2,4),hydrogen(j2),helium4(j2),helium3(j2),iflip)
        call nabla(X(J3,1),X(J3,2),X(J3,3),X(J3,4),zMhf3,               &
     &     TPNAB3,ADNAB,RADNAB,J3,IC)
        iconvect(j3) = ic 
        QNh1 = hydrogen(j2)*RHO2/HMASS
        QNh2 =deuterium(j2)*RHO2/2.d0/HMASS
        QNhe3 = helium3(j2)*RHO2/3.d0/HMASS
        CALL NUCRAT(RHO2,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EP,              &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
        if(K2.gt.0) then
          dNAB2P = dNAB3P
          dNAB2R = dNAB3R
          dNAB2L = dNAB3L
          dNAB2T = dNAB3T
          TPNAB2P=TPNAB3P
          TPNAB2R=TPNAB3R
          TPNAB2L=TPNAB3L
          TPNAB2T=TPNAB3T
          PN=X(J2,1)*(1.d0+XINC)
          PM=X(J3,1)*(1.d0+XINC)
          call invstate(27,PN,x(j2,4),hydrogen(j2),helium4(j2),Crad,    &
     &                  RHP2,cPP,alphaP,beta,deltaP)
          call nabla(PM,X(J3,2),X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3P,ADNAB,RADNAB,J3,IC)
          RM=X(J3,2)*(1.d0+XINC)
          call nabla(X(J3,1),RM,X(J3,3),X(J3,4),zMhf3,                  &
     &     TPNAB3R,ADNAB,RADNAB,J3,IC)
          XLM=X(J3,3)*(1.d0+XINC)
          call nabla(X(J3,1),X(J3,2),XLM,X(J3,4),zMhf3,                 &
     &     TPNAB3L,ADNAB,RADNAB,J3,IC)
          TN=X(J2,4)*(1.d0-XINC)
          TM=X(J3,4)*(1.d0-XINC)
          call invstate(28,x(j2,1),TN,hydrogen(j2),helium4(j2),Crad,    &
     &                  RHT2,cPT,alphaT,beta,deltaT)
          CALL NUCRAT(RHO2,TN,QNh1,QNh2,QNhe3,DTIME,EPT,                &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          call nabla(X(J3,1),X(J3,2),X(J3,3),TM,zMhf3,                  &
     &     TPNAB3T,ADNAB,RADNAB,J3,IC)
          RHD=RHO2*(1.d0+XINC)
          QNh1 = hydrogen(j2)*rhD /HMASS
          QNh2 =deuterium(j2)*rhD /2.d0/HMASS
          QNhe3 = helium3(j2)*rhD /3.d0/HMASS
          CALL NUCRAT(RHD,X(J2,4),QNh1,QNh2,QNhe3,DTIME,EPD,            &
     &      Rpp,Rhe3,Rcno,R3a,R12a,Xpd,Xpc,Xpn,Xpo,                     &
     &      helium4(j2),carbon(j2),trogen(j2),oxygen(j2),               &
     &      xpp23,iflip,-j2)
          dRHdP2=(RHP2-RHO2)/(XINC*X(J2,1))
          dRHdT2=(RHT2-RHO2)/(-XINC*X(J2,4))
          dEPdR=(EPD-EP)/(XINC*RHO2)
          dEPdT=(EPT-EP)/(-XINC*X(J2,4)) + dEPdR*dRHdT2
          dEPdP=dEPdR*dRHdP2
          dcPdP=(cPP-cP)/(XINC*X(J2,1))
          dcPdT=(cPT-cP)/(-XINC*X(J2,4))
          ddeldP=(deltaP-delta)/(XINC*X(J2,1))
          ddeldT=(deltaT-delta)/(-XINC*X(J2,4))
          dNAB3P=(TPNAB3P-TPNAB3)/(XINC*X(J3,1))
          dNAB3R=(TPNAB3R-TPNAB3)/(XINC*X(J3,2))
          dNAB3L=(TPNAB3L-TPNAB3)/(XINC*X(J3,3))
          dNAB3T=(TPNAB3T-TPNAB3)/(-XINC*X(J3,4))
          HD(1,1) =-1.d0
          HD(1,2) =-4.d0*GMdM/X(J2,2)**5
          HE(1,1) = 1.d0
          HC(2,2) =-3.d0*X(J1,2)**2
          HD(2,1) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdP2
          HD(2,2) = 3.d0*X(J2,2)**2
          HD(2,4) = dM(J2)*0.75d0/(PI*RHO2**2)*dRHdT2
          HC(3,3) =-1.d0
          HD(3,1) =-dM(J2)*( dEPdP                                      &
     &             -dcPdP * (X(J2,4)-VX(J2,4))*brem                     &
     &             +delta/RHO2 *brem                                    &
     &             +ddeldP/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdP2 * (X(J2,1)-VX(J2,1))*brem )
          HD(3,3) = 1.d0
          HD(3,4) =-dM(J2)*( dEPdT                                      &
     &             -cP * brem                                           &
     &             -dcPdT * (X(J2,4)-VX(J2,4))*brem                     &
     &             +ddeldT/RHO2 * (X(J2,1)-VX(J2,1))*brem               &
     &             -delta/RHO2**2 * dRHdT2 * (X(J2,1)-VX(J2,1))*brem )
          HD(4,1) = 0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB2P
          HD(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2R
          HD(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB2L
          HD(4,4) =-1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB2T
          HE(4,1) =-0.5d0*(TPNAB3+TPNAB2)                               &
     &                   - (X(J3,1)-X(J2,1))*0.5d0*dNAB3P
          HE(4,2) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3R
          HE(4,3) =      - (X(J3,1)-X(J2,1))*0.5d0*dNAB3L
          HE(4,4) = 1.d0 - (X(J3,1)-X(J2,1))*0.5d0*dNAB3T
        endif
        G(1) = X(J3,1)-X(J2,1) + GMdM/X(J2,2)**4
        G(2) = X(J2,2)**3 - X(J1,2)**3 - dM(J2)*0.75d0/(PI*RHO2)
        G(3) = X(J2,3) - X(J1,3) - dM(J2)*( EP                          &
     &       - cP * (X(J2,4)-VX(J2,4))*brem                             &
     &       + delta/RHO2 * (X(J2,1)-VX(J2,1))*brem )
        G(4) = (X(J3,4)-X(J2,4))                                        &
     &         -(X(J3,1)-X(J2,1))*0.5d0*(TPNAB3+TPNAB2)
C       write(6,*) j2, ' Gs ', g(1),g(2),g(3),g(4) 
        return
      endif
      if (J2.eq.N) then
        if(iter.ge.abs(itmax/2).or.X(j2,3).lt.0.) then
          call atmos(Tatm,RHOatm,Ratm,Patm,VX(J2,2),VX(J2,3),0)
          if(K2.gt.0) then
            HD(1,1) = 1.d0
            HD(2,2) = 1.d0
            HC(2,2) =-1.d0
            HD(3,3) = 1.d0
            HC(3,3) =-1.d0
            HD(4,4) = 1.d0
          endif
        else
          call atmos(Tatm,RHOatm,Ratm,Patm,X(J2,2),X(J2,3),0)
          if(K2.gt.0) then
            RN=X(J2,2)*(1.d0+XINC )
            call atmos(TatR,RHOatR,RatR,PatR,RN,X(J2,3),0)
            xLN=X(J2,3)*(1.d0+XINC) 
            call atmos(TatL,RHOatL,RatL,PatL,X(J2,2),xLN,0)
            dPadr = (PatR - Patm)/(XINC * X(J2,2))
            dPadL = (PatL - Patm)/(XINC * X(J2,3))
            dTadr = (TatR - Tatm)/(XINC * X(J2,2))
            dTadL = (TatL - Tatm)/(XINC * X(J2,3))
            dRadr = (RHOatR - RHOatm)/(XINC * X(J2,2))
            dRadL = (RHOatL - RHOatm)/(XINC * X(J2,3))
            datrdr = (RatR - Ratm)/(XINC * X(J2,2))
            datrdL = (RatL - Ratm)/(XINC * X(J2,3))
c           write(6,*) ' atmos derivs PTR ', dpadr,dpadl,dtadr,dtadl,   &
c    &                 dradr,dradl 
            HD(1,1) = 1.d0
            HD(1,2) = -dPadr
            HD(1,3) = -dPadL
            HD(2,2) = 1.d0 - dredgedr - datrdR
            HD(2,3) =      - dredgedl - datrdL
            HC(2,2) =-1.d0
            HD(3,3) = 1.d0
            HC(3,3) =-1.d0
            HD(4,4) = 1.d0
            HD(4,2) = -dTadr
            HD(4,3) = -dTadL
          endif
        endif
        G(1) = X(J2,1)-Patm
        G(2) = X(J2,2)-X(J1,2) - Ratm
        G(3) = X(J2,3)-X(J1,3)
        G(4) = X(J2,4)-Tatm
c       write(6,*) 'atmos Gs ', g(1),g(2),g(3),g(4) 
        return
      endif 
      return
      end
c
      subroutine nabla(PRES,RAD,XLUM,TEMP,XMASS,TPNAB,ADNAB,RADNAB,J,IC)
C
C ** CHECK MIXING LENGTH GRADIENT AND CONVECTIVE VELOCITIES. **
C ** KIPPENHAHN,WEIGERT,HOFMEISTER... AND TALBOT AND SMITH
C
      include 'parm.h'
      include 'var.h'
c     DATA anabcon/3.945322d+09/
      DATA GRAV/6.6704D-8/
      DATA PI  /3.14159265359d0/
      DATA SIG,RG/5.67051d-5,8.31451d7/
      data iwrite/0/
      data Trans,TRwidth/4.d5,10.d0/
      save

      anabcon=3.d0/(GRAV*PI*64.d0*SIG)
      K=abs(J)
        call invstate(31,PRES,TEMP,hydrogen(K),helium4(K),Crad,         &
     &     RHO,cP,alpha,beta,delta)
        RHOlog=log10(RHO)
        Tlog  =log10(TEMP)
        call opacity(K,Tlog,RHOlog,bkap) 
        xkap=10.d0**bkap
        radnab=anabcon/xMASS*xkap*xLUM*PRES/TEMP**4
        adnab=delta*PRES/(cP*TEMP*RHO)
       if(radnab.le.adnab) then
         IC=1
         TPNAB=TEMP/PRES*radnab
       else
         IC=0
         TPNAB=TEMP/PRES*adnab
         FACTR=max(0.d0,(vx(k,4)/Trans-1.d0)*TRwidth)
        if(FACTR.gt.1.d0) return
         WTM=Rg*RHO*TEMP/PRES
         GR = GRAV*xmass/RAD**2
         BET=1.d0
         HP=Rg*TEMP/GR/WTM/BET
         DIST=RLH*HP
C      CONVECTIVE GRADIENT
         A1=12.d0*SIG*TEMP**3/CP/RHO/RHO/DIST/XKAP
         YOU=A1*SQRT(8.d0*HP/GR/DIST/DIST/DELTA)
         E3=0.5048011d0*YOU*YOU
         A2=YOU*(4.d0/9.d0*(radnab-adnab) + 0.2373622d0*YOU*YOU)
         ddd=A2*A2+E3*E3*E3
         Wd=(A2+sqrt(ddd))**(1.d0/3.d0)
         C1=(Wd*Wd+0.7037037*Wd*YOU-E3)**2
         DCONV=max(0.d0,C1/Wd/Wd-YOU*YOU)+adnab
         if (dconv.lt. 0.d0)then
           dconv=radnab
         end if
         TPNAB=FACTR*TPNAB + (1.d0-FACTR)*TEMP/PRES*dconv
       end if
      return
      end
      subroutine GID(J2)
C+
C  The subroutine GID calculates the numerical derivatives to the
C  equations given in GI and compares them to the analytical ones.
C
C  Author:  Harold W. Yorke (JPL)
C  Date:    29-Jan-2003
C-
      include 'parm.h'
      include 'var.h'
      dimension gk(MH),gg(MH)
      data DINC/-0.0001/
c
      save
c
      j1=j2-1
      j3=j2+1
      write(6,221) '***********  GI(',j2,')=',(g(i),i=1,MH)
 221  format(a,i4,a,1p,4E12.4)
 222  format(i5,1p,4E10.2)
 223  format(5x,1p,4E10.2)
 224  format(5x,1p,4E10.2,' ********')
      do i=1,MH
        gk(i)=g(i)
      enddo
      call gi(-j2)
      write(6,221) ' numerical   GI(',j2,')=',(g(i),i=1,MH)
      if(j2.gt.1) then
      do k=1,MH
        xkeep=x(j1,k)
        x(j1,k) = x(j1,k)*(1.+DINC)
        call gi(-j2)
        x(j1,k) = xkeep
        do i=1,MH
          gg(i)=(g(i)-gk(i))/(DINC*xkeep+1.E-37)
        enddo
        write(6,222) j1,(gg(i),i=1,MH)
        write(6,222) k,(HC(i,k),i=1,MH)
        dddx=0.
        do i=1,MH
          xxxx = gg(i)
          gg(i)=.5*(HC(i,k)-gg(i))/(1.E-37+HC(i,k)+gg(i))
          HC(i,k)=xxxx
          dddx=dddx+abs(gg(i))
        enddo
        if(dddx.lt.0.2) then
          write(6,223) (gg(i),i=1,MH)
        else
          write(6,224) (gg(i),i=1,MH)
        endif
      enddo
      endif
      do k=1,MH
        xkeep = x(j2,k)
        x(j2,k) = x(j2,k)*(1.+DINC)
        call gi(-j2)
        x(j2,k) = xkeep
        do i=1,MH
          gg(i)=(g(i)-gk(i))/(DINC*xkeep+1.E-37)
        enddo
        write(6,222) j2,(gg(i),i=1,MH)
        write(6,222) k,(HD(i,k),i=1,MH)
        dddx=0.
        do i=1,MH
          xxxx = gg(i)
          gg(i)=.5*(HD(i,k)-gg(i))/(1.E-37+HD(i,k)+gg(i))
          HD(i,k)=xxxx
          dddx=dddx+abs(gg(i))
        enddo
        if(dddx.lt.0.2) then
          write(6,223) (gg(i),i=1,MH)
        else
          write(6,224) (gg(i),i=1,MH)
        endif
      enddo
      if(j2.lt.N) then
      do k=1,MH
        xkeep = x(j3,k)
        x(j3,k) = x(j3,k)*(1.+DINC)
        call gi(-j2)
        x(j3,k) = xkeep
        do i=1,MH
          gg(i)=(g(i)-gk(i))/(DINC*xkeep+1.E-37)
        enddo
        write(6,222) j3,(gg(i),i=1,MH)
        write(6,222) k,(HE(i,k),i=1,MH)
        dddx=0.
        do i=1,MH
          xxxx = gg(i)
          gg(i)=.5*(HE(i,k)-gg(i))/(1.E-37+HE(i,k)+gg(i))
          HE(i,k)=xxxx
          dddx=dddx+abs(gg(i))
        enddo
        if(dddx.lt.0.2) then
          write(6,223) (gg(i),i=1,MH)
        else
          write(6,224) (gg(i),i=1,MH)
        endif
      enddo
      endif
      call gi(j2)
      do i=1,MH
        g(i)=gk(i)
      enddo
      return
      end
