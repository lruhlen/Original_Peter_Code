      SUBROUTINE ADDSUB
C+
C  THE SUBROUTINE ADDSUB DETERMINES WHETHER A GRID POINT SHOULD BE ADDED
C  OR REMOVED.
C
C  Author:  H.W. YORKE   23-SEP-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C  
C-
      include 'parm.h'
      include 'var.h'
      parameter(Mk=17)
      character*1 why(3)
      character*3 blank,why3
      equivalence (why,why3)
      data blank/'   '/
      data Nadd,Nsub/0,0/
      save
c
c CHECK WHETHER A POINT SHOULD BE ADDED TO THE OUTER ZONE
c
      why3=blank
      Jadd=0
      Jsub=0
      if(x(N,4).gt.Atmx) then
        do jjj=1,2
          Jadd=N
          why(1)='A'
          call add(Jadd,why)
          Nadd=Nadd+1
        enddo
      endif
c
c CHECK WHETHER A POINT SHOULD BE REMOVED FROM THE OUTER ZONE
c
      if(x(N,4).lt.Atmn .and. (dM(N-1)+dM(N))/zM(N-1) .lt. dzmax) then
        Jsub=N-1
        call sub(Jsub)
        Nsub=Nsub+1
      endif
c
c CHECK WHETHER MORE POINTS SHOULD BE ADDED CLOSE TO THE EXTERIOR
c
      if(Jadd+Jsub .eq. 0) then
        do j=N-5,N-1
          if(dM(j)/zM(j) .gt. dzmax) then
            Jadd = J
            why(1)='Z'
            call add(Jadd,why)
            Nadd=Nadd+1
            goto 10
          endif
        enddo
      endif
 10   continue
c
c CHECK WHETHER GRID POINTS SHOULD BE ADDED OR REMOVED
c
      dPL =1.E37
      Slumax=0.
      do j=1,N
        Slumax=max(Slumax,abs(x(j,3)))
      enddo
      Jadd=0
      Jsub=0
      do j=2,99999
        dP   =abs((x(j,1)-x(j-1,1))/(x(j,1)+x(j-1,1)))
        dLum =abs((x(j,3)-x(j-1,3))/Slumax)
        dX   =abs(hydrogen(j)-hydrogen(j-1))
     *       +abs(helium4 (j)-helium4 (j-1))
        if(J.le.N-2) then
          dP2  =abs((x(j+1,1)-x(j,1))/(x(j,1)+x(j+1,1)))
          dLum2=abs((x(j+1,3)-x(j,3))/Slumax)
          dX2  =abs(hydrogen(j)-hydrogen(j+1))
     *         +abs(helium4 (j)-helium4 (j+1))
          dP3  =abs((x(j+1,1)-x(j+2,1))/(x(j+2,1)+x(j+1,1)))
          dLum3=abs((x(j+1,3)-x(j+2,3))/Slumax)
          dX3  =abs(hydrogen(j+2)-hydrogen(j+1))
     *         +abs(helium4 (j+2)-helium4 (j+1))
          if((J-Jsub.gt.2 .and. J-Jadd.gt.2 .and.
     *      dLum.lt.dLmn .and. dLum2.lt.dLmn .and. dLum3.lt.dLmn .and.
     *      dP.lt.dPmn .and. dP2.lt.dPmn .and. dP3.lt.dPmn .and.
     *      dX.lt.dXmn .and. dX2.lt.dXmn .and. dX3.lt.dXmn .and.
     *      (dM(j+1)+dM(j))/zM(N) .lt. dzmax)
     *      .or. (dM(j+1)+dM(j))/zM(j) .lt. dzmin)
     *    then
            Jsub=j
            call sub(Jsub)
            Nsub=Nsub+1
          endif
        endif
        iwhy=0
        why3=blank
        if(dLum.gt.dLmx) then
          iwhy=iwhy+1
          why(iwhy)='L'
        endif
        if(dP.gt.dPmx) then
          iwhy=iwhy+1
          why(iwhy)='P'
        endif
        if(dX.gt.dXmx) then
          iwhy=iwhy+1
          why(iwhy)='X'
        endif
        if(J-Jadd.gt.2 .and. iwhy.ne.0 .and. N.lt.MJ .and.
     *     dM(J)/zM(J) .gt. 2.*dzmin .and. J-Jsub.gt.1) then
          Jadd = J
          call add(Jadd,why)
          Nadd=Nadd+1
        endif
        if(j.ge.N-1) goto 20
      enddo
 20   continue
      if(Nadd+Nsub.gt.0 .and. mod(MODEL,JRIT).eq.0 .and. JRIT.gt.1) then
        write(6,201) Nsub,Nadd,N
 201    format(' ===> A total of ',I4,' points were subtracted and ',I4,
     *         ' points were added. N=',I4)
        Nadd=0
        Nsub=0
      endif
      return
      end
      SUBROUTINE ADD(J,why)
C+
C  THE SUBROUTINE ADD ADDS A GRID POINT AT J
C
C  Author:  H.W. YORKE   30-AUG-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      include 'var.h'
      character*1 why(3)
      data PI/3.14159265359d0/
      save
C
      if(N.lt.MJ .and. J.eq.N) then
        dM(J)=dM(J)*0.5d0
        dM(J+1)=dM(J)
        zM(J+1)=zM(J)
        zM(J) = zM(J)-dM(J)
        hydrogen(J+1) = hydrogen(J)
        helium3 (J+1) = helium3 (J)
        helium4 (J+1) = helium4 (J)
        deuterium(J+1)=deuterium(J)
        carbon  (J+1) = carbon  (J)
        trogen  (J+1) = trogen  (J)
        oxygen  (J+1) = oxygen  (J)
        vhydrogen(J+1) = vhydrogen(J)
        vhelium3 (J+1) = vhelium3 (J)
        vhelium4 (J+1) = vhelium4 (J)
        vdeuterium(J+1)=vdeuterium(J)
        vcarbon  (J+1) = vcarbon  (J)
        vtrogen  (J+1) = vtrogen  (J)
        voxygen  (J+1) = voxygen  (J)
C  Interpolate/extrapolate for new values of pressure at J=N,N+1
        w1 = dM(J)/(dM(J-1)+dM(J)+dM(J+1))
        dP = (x(J,1)-x(J-1,1))*w1
        x (J+1,1) = x (J,1) + dP
        x (J  ,1) = x (J,1) - dP
        dP = (vx(J,1)-vx(J-1,1))*w1
        vx(J+1,1) = vx(J,1) + dP
        vx(J  ,1) = vx(J,1) - dP
C  Use average density to estimate radius of new gridpoint at J=N
        x(J+1,2) = x(J,2)
        RHOPI  = (x(J,2)**3-x(J-1,2)**3)/(dM(J)+dM(J+1))
        x(J,2) = (x(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
        vx(J+1,2) = vx(J,2)
        RHOPI  = (vx(J,2)**3-vx(J-1,2)**3)/(dM(J)+dM(J+1))
        vx(J,2)= (vx(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
C  Interpolate for new values of luminosity at J=N
        x (J+1,3) = x(J,3)
        vx(J+1,3) = vx(J,3)
        x (J,3) = 0.5d0*(x(J-1,3)+x(J,3))
        vx(J,3) = 0.5d0*(vx(J-1,3)+vx(J,3))
C  Interpolate/extrapolate for new values of temperature at J=N,N+1
        dT = (x(J,4)-x(J-1,4))*w1
        x (J+1,4) = x (J,4) + dT
        x (J  ,4) = x (J,4) - dT
        dT = (vx(J,1)-vx(J-1,1))*w1
        vx(J+1,4) = vx(J,4) + dT
        vx(J  ,4) = vx(J,4) - dT
      endif
      if(N.lt.MJ .and. J.lt.N) then
        do i=N,J,-1
          dM      (i+1) = dM      (i)
          zM      (i+1) = zM      (i)
          hydrogen(i+1) = hydrogen(i)
          helium3 (i+1) = helium3 (i)
          helium4 (i+1) = helium4 (i)
          deuterium(i+1)=deuterium(i)
          carbon  (i+1) = carbon  (i)
          trogen  (i+1) = trogen  (i)
          oxygen  (i+1) = oxygen  (i)
          vhydrogen(i+1) = vhydrogen(i)
          vhelium3 (i+1) = vhelium3 (i)
          vhelium4 (i+1) = vhelium4 (i)
          vdeuterium(i+1)=vdeuterium(i)
          vcarbon  (i+1) = vcarbon  (i)
          vtrogen  (i+1) = vtrogen  (i)
          voxygen  (i+1) = voxygen  (i)
          do k=1,MH
            vx(i+1,k)=vx(i,k)
            x (i+1,k)=x (i,k)
          enddo
        enddo
        dM(J)=dM(J)*0.5d0
        dM(J+1)=dM(J)
        zM(J) = zM(J)-dM(J)
        if(J.gt.1) then
C  Interpolate for new values of pressure
          w1 = dM(J)/(4.d0*dM(J)+dM(J-1)+dM(J+2))
          dP = (x(J+2,1)-x(J-1,1))*w1
          x (J  ,1) = x (J  ,1) - dP
          x (J+1,1) = x (J+1,1) + dP
          dP = (vx(J+2,1)-vx(J-1,1))*w1
          vx(J  ,1) = vx(J  ,1) - dP
          vx(J+1,1) = vx(J+1,1) + dP
C  Use average density to estimate radius of new gridpoint
          RHOPI  = (x(J+2,2)**3-x(J-1,2)**3)/(dM(J)+dM(J+1))
          x(J,2) = (x(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
          RHOPI  = (vx(J+2,2)**3-vx(J-1,2)**3)/(dM(J)+dM(J+1))
          vx(J,2)= (vx(J-1,2)**3 + dM(J)*RHOPI)**(1./3.)
C  Interpolate for new values of luminosity
          x (J,3) = 0.5d0*(x(J-1,3)+x(J+1,3))
          vx(J,3) = 0.5d0*(vx(J-1,3)+vx(J+1,3))
C  Interpolate for new values of temperature
          dT = (x(J+2,4)-x(J-1,4))*w1
          x (J  ,4) = x (J  ,4) - dT
          x (J+1,4) = x (J+1,4) + dT
          dT = (vx(J+2,1)-vx(J-1,1))*w1
          vx(J  ,4) = vx(J  ,4) - dT
          vx(J+1,4) = vx(J+1,4) + dT
        else
C  Interpolate for new values of pressure at J=1
          w1 = dM(J)/(dM(J)+dM(J+1)+dM(J+2))
          dP = (x(J+2,1)-x(J+1,1))*w1
          x (J  ,1) = x (J  ,1) - dP
          x (J+1,1) = x (J+1,1) + dP
          dP = (vx(J+2,1)-vx(J+1,1))*w1
          vx(J  ,1) = vx(J  ,1) - dP
          vx(J+1,1) = vx(J+1,1) + dP
C  Use average density to estimate radius of new gridpoint at J=1
          RHOPI  = x(J+2,2)**3/(dM(J)+dM(J+1))
          x(J,2) = (dM(J)*RHOPI)**(1./3.)
          RHOPI  = vx(J+2,2)**3/(dM(J)+dM(J+1))
          vx(J,2)= (dM(J)*RHOPI)**(1./3.)
C  Interpolate for new values of luminosity at J=1
          x (J,3) = 0.5d0*x(J+1,3)
          vx(J,3) = 0.5d0*vx(J+1,3)
C  Interpolate for new values of temperature at J=1
          dT = (x(J+2,4)-x(J+1,4))*w1
          x (J  ,4) = x (J  ,4) - dT
          x (J+1,4) = x (J+1,4) + dT
          dT = (vx(J+2,1)-vx(J+1,1))*w1
          vx(J  ,4) = vx(J  ,4) - dT
          vx(J+1,4) = vx(J+1,4) + dT
        endif
      endif
      if(N.lt.MJ) then
        N=N+1
        if(mod(MODEL,JRIT).eq.0 .and. JRIT.eq.1)
     *    write(6,201) J,zM(J),dM(J+1),why
 201      format(' ===> New grid point added at J=',i4,'. M=',1p,E13.6,
     *           ' dM=',E11.4,1x,3a1)
      endif
      return
      end
      SUBROUTINE SUB(J)
C+
C  THE SUBROUTINE SUB REMOVES A GRID POINT AT J
C
C  Author:  H.W. YORKE   23-AUG-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      include 'var.h'
      save
C
      if(J.lt.N .and. J.ge.1) then
        dMJ=dM(J)+dM(J+1)
        w2=dM(J)/dMJ
        w1=dM(J+1)/dMJ
        hydrogen(J) = w1*hydrogen(J) + w2*hydrogen(J+1)
        helium3 (J) = w1*helium3 (J) + w2*helium3 (J+1)
        helium4 (J) = w1*helium4 (J) + w2*helium4 (J+1)
        deuterium(J)= w1*deuterium(J)+ w2*deuterium(J+1)
        carbon  (J) = w1*carbon  (J) + w2*carbon  (J+1)
        trogen  (J) = w1*trogen  (J) + w2*trogen  (J+1)
        oxygen  (J) = w1*oxygen  (J) + w2*oxygen  (J+1)
        vhydrogen(J) = w1*vhydrogen(J) + w2*vhydrogen(J+1)
        vhelium3 (J) = w1*vhelium3 (J) + w2*vhelium3 (J+1)
        vhelium4 (J) = w1*vhelium4 (J) + w2*vhelium4 (J+1)
        vdeuterium(J)= w1*vdeuterium(J)+ w2*vdeuterium(J+1)
        vcarbon  (J) = w1*vcarbon  (J) + w2*vcarbon  (J+1)
        vtrogen  (J) = w1*vtrogen  (J) + w2*vtrogen  (J+1)
        voxygen  (J) = w1*voxygen  (J) + w2*voxygen  (J+1)
        x (J,1) = w1*x(J,1) + w2*x(J+1,1)
        x (J,2) = x(J+1,2)
        x (J,3) = x(J+1,3)
        x (J,4) = w1*x(J,4) + w2*x(J+1,4)
        vx(J,1) = w1*vx(J,1) + w2*vx(J+1,1)
        vx(J,2) = x(J+1,2)
        vx(J,3) = x(J+1,3)
        vx(J,4) = w1*vx(J,4) + w2*vx(J+1,4)
        dM(J) = dMJ
        zM(J)  = zM(J+1)
        if(J+1.lt.N-1) then
          do i=J+1,N-1
            zM      (i) = zM      (i+1)
            dM      (i) = dM      (i+1)
            hydrogen(i) = hydrogen(i+1)
            helium3 (i) = helium3 (i+1)
            helium4 (i) = helium4 (i+1)
            deuterium(i)=deuterium(i+1)
            carbon  (i) = carbon  (i+1)
            trogen  (i) = trogen  (i+1)
            oxygen  (i) = oxygen  (i+1)
            vhydrogen(i) = vhydrogen(i+1)
            vhelium3 (i) = vhelium3 (i+1)
            vhelium4 (i) = vhelium4 (i+1)
            vdeuterium(i)=vdeuterium(i+1)
            vcarbon  (i) = vcarbon  (i+1)
            vtrogen  (i) = vtrogen  (i+1)
            voxygen  (i) = voxygen  (i+1)
            do k=1,MH
              vx(i,k)=vx(i+1,k)
              x (i,k)=x (i+1,k)
            enddo
          enddo
        endif
        N=N-1
      endif
      if(JRIT.eq.1) then
        if(J.eq.1) then
          if(mod(MODEL,JRIT).eq.0) write(6,201) J,zM(J)
        else
          if(mod(MODEL,JRIT).eq.0) write(6,201) J,zM(J-1),zM(J)
        endif
      endif
 201  format(' ===> Grid point removed at J=',i4,'. M=',1p,E13.6,       &
     &     ' dM=',E11.4)
      return
      end
      SUBROUTINE GRIDMOV
C+
C  THE SUBROUTINE GRIDMOV MODIFIES THE MASS DISTRIBUTION OF GRID
C
C  Author:  H.W. YORKE   18-DEC-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      include 'var.h'
      parameter(Mk=15)
      dimension WT(Mk),Zold(Mk),Fold(Mk+2),z(Mk)
      save

      Znorm = zM(N)-zM(N-Mk+1)
      Zold(1)=0.d0
      do i=2,Mk
        j=N-Mk+i
        Zold(i)=Zold(i-1)+dM(j)/Znorm
      enddo
C
C  Try to keep the bottom of the atmosphere at Tmean (Tmean is defined as
C  the average of the allowed maximum and minimum atmasphere temperatures).
      Tmean=.5*(Atmx+Atmn)
      dZ=X(N,4)/Tmean-1.d0
      if(dZ.gt.0.d0) then
        dZ=min(.5d0,dZ)*dZmax
      else
        dZ=max(-.5d0,dZ)*dZmax
      endif
      dZ=dZ*dZdt
      WX = min(.5d0,dZdt*0.05d0)
      do i=2,Mk-1
        zev  = (1.d0+dZ)*Zold(Mk-1)*float(i-1)/float(Mk-2)
        ddz  = min(0.5d0*(Zold(i+1)-Zold(i)),zev-Zold(i))
        ddz  = max(ddz,0.5d0*(Zold(i-1)-Zold(i)))
        z(i) = Zold(i) + WX*ddz
      enddo
      z(1) = 0.d0
      z(Mk)= 1.d0
      do i=2,Mk-1
        if(z(i).lt.Zold(i)) then
          WT(i)=(z(i)-Zold(i))/(Zold(i)-Zold(i-1))
        else
          WT(i)=(z(i)-Zold(i))/(Zold(i+1)-Zold(i))
        endif
      enddo
      WT(1)=0.d0
      WT(Mk)=0.d0
C  Interpolate quantities defined on mass grid (R and L)
      call interpol(WT,x (1,2),Fold,N,Mk)
      call interpol(WT,vx(1,2),Fold,N,Mk)
      call interpol(WT,x (1,3),Fold,N,Mk)
      call interpol(WT,vx(1,3),Fold,N,Mk)
C  Redefine mass grid 
      zMass = zM(N)
      do i=2,Mk
        j=N-Mk+i
        dM(j) = Znorm*(z(i)-z(i-1))
        zM(j)  = zM(j-1) + dM(j)
      enddo
      zM(N) = zMass
      if(mod(MODEL,JRIT).eq.0)
     *     write(6,201) N,(z(i)-Zold(i),i=Mk-3,Mk-1),(x(j,4),j=N-2,N)
 201  format(' GRIDMOV:  N=',i4,1p,6E11.3)
C  Redefine interpolation grid to middle of mass grid
      do i=2,Mk
        z(i)=0.5d0*(z(i-1)+z(i))
        Zold(i)=0.5d0*(Zold(i-1)+Zold(i))
      enddo
      Zold(1) = -0.5d0*dM(N-Mk+1)/Znorm
      do i=2,Mk-1
        if(z(i).lt.Zold(i)) then
          WT(i)=(z(i)-Zold(i))/(Zold(i)-Zold(i-1))
        else
          WT(i)=(z(i)-Zold(i))/(Zold(i+1)-Zold(i))
        endif
        WT(Mk)=min(0.d0,(z(Mk)-Zold(Mk))/(Zold(Mk)-Zold(Mk-1)))
      enddo
C  Interpolate quantities defined between mass gridpoints
      call interpol(WT,x (1,1),Fold,N,Mk)
      call interpol(WT,vx(1,1),Fold,N,Mk)
      call interpol(WT,x (1,4),Fold,N,Mk)
      call interpol(WT,vx(1,4),Fold,N,Mk)
      call interpol(WT,hydrogen,Fold,N,Mk)
      call interpol(WT,helium3 ,Fold,N,Mk)
      call interpol(WT,helium4 ,Fold,N,Mk)
      call interpol(WT,deuterium,Fold,N,Mk)
      call interpol(WT,carbon  ,Fold,N,Mk)
      call interpol(WT,trogen  ,Fold,N,Mk)
      call interpol(WT,oxygen  ,Fold,N,Mk)
      call interpol(WT,vhydrogen,Fold,N,Mk)
      call interpol(WT,vhelium3 ,Fold,N,Mk)
      call interpol(WT,vhelium4 ,Fold,N,Mk)
      call interpol(WT,vdeuterium,Fold,N,Mk)
      call interpol(WT,vcarbon  ,Fold,N,Mk)
      call interpol(WT,vtrogen  ,Fold,N,Mk)
      call interpol(WT,voxygen  ,Fold,N,Mk)
      return
      end
      subroutine interpol(WT,F,Fold,N,Mk)
      include 'parm.h'
      dimension WT(Mk),F(N),Fold(Mk+2)
      save

      do i=1,Mk
        j=N-Mk+i
        Fold(i)=F(j)
      enddo
      do i=2,Mk-1
        if(WT(i).gt.0.d0) then
          W3=WT(i)
          W2=1.d0-W3
          W1=0.d0
        else
          W1=-WT(i)
          W2=1.d0-W1
          W3=0.d0
        endif
        j=N-Mk+i
        F(j) = W1*Fold(i-1)+W2*Fold(i)+W3*Fold(i+1)
      enddo
      W1=-WT(Mk)
      W2=1.d0-W1
      F(N) = W1*Fold(Mk-1) + W2*Fold(Mk)
      return
      end
      SUBROUTINE MASSFLUX(time,Zflux)
C+
C  THE SUBROUTINE MASSFLUX SPECIFIES THE CURRENT MASS FLUX ONTO
C  THE STAR (OR MASS LOSS FROM THE STAR IF ZFLUX IS NEGATIVE)
C
C  Author:  H.W. YORKE   18-DEC-02               (JPL / CALTECH)
C  Version: 12-May-2006
C  History: (see summary in CHANGES)
C-
      include 'parm.h'
      parameter (MF=20)
      dimension tim(MF),flx(MF)
      data ifirst/0/
      save

      if(ifirst.eq.0) then
        ifirst=1
        do i=1,MF
          tim(i)=-1.
          call comrd
          READ (5,*,err=10,end=10) TIM(i),FLX(i)
          if(tim(i).lt.0.d0) goto 10
        enddo
        i=MF
 10     NF=i-1
      endif
      do i=1,NF
        if(time.lt.tim(i)) goto 20
      enddo
      i=NF+1
 20   continue
      if(i.eq.1) then
        Zflux=flx(1)
      else
        if(i.gt.NF) then
          Zflux=flx(NF)
        else
          dt=(time-tim(i-1))/(tim(i)-tim(i-1))
          Zflux=flx(i-1) + dt * (flx(i)-flx(i-1))
        endif
      endif
      return
      end
      subroutine addmass
      include 'parm.h'
      include 'var.h'
      save
c
      call massflux(time,Zflux)
      dM(N) = dM(N) + Zflux*dtime
      if(dM(N).le.0.d0) stop 'addmass'
c     do j=1,N
c       dM(j) = dM(j)*(1.+1.d-6)
c       zM(j) = zM(j)*(1.+1.d-6)
c     enddo
      zM(N)  = zM(N-1) + dM(N)
      Zmass  = zM(N)
      return
      end
