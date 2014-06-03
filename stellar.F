      PROGRAM STELLAR
      include 'parm.h'
      include 'var.h'

      save
c
      call start
      imod=0
      konv=0
      mode=0
      do i=1,100000
        imod=imod+1
        TIME=TIME+DTIME
        MODEL=MODEL+1
C
C  Store old values of time dependent variables in arrays with names
C  beginning with a "v".
C
        do j=1,N
          do k=1,NG
            vx(j,k)=x(j,k)
          enddo
          vhydrogen(j)=hydrogen(j)
          vhelium4(j)=helium4(j)
          vhelium3(j)=helium3(j)
          vdeuterium(j)=deuterium(j)
          vcarbon(j)=carbon(j)
          vtrogen(j)=trogen(j)
          voxygen(j)=oxygen(j)
        enddo
C
C  End of storing old values of time dependent variables
C
        konv=max(0,konv-1)
        call massflux(TIME,Zflux)
        if(konv.gt.0. or. MODEL.eq.1) then
          mode=999
        else
          mode=0
          call gridmov
        endif
        dTAUmx= taufactor*taujk
c       call elatim('NOPRINT')
        call henyey(mode)
        if(mode.ne.0) then
          konv=konv+3
C
C  NO CONVERGENCE IN HENYEY.  TRY AGAIN
C
          TIME=TIME-DTIME
          MODEL=MODEL-1
          imod=imod-1
          if(abs(DTIME/DTMIN-1.).lt.0.01) stop 'NO CONVERGENCE'
          DTIME=DTIME/float(konv)
          DTIME=max(DTIME,DTMIN)
C
C  Reset time dependent variables with old values.
C
          do j=1,N
            do k=1,NG
              x(j,k)=vx(j,k)
            enddo
            hydrogen(j)=vhydrogen(j)
            helium4(j)=vhelium4(j)
            helium3(j)=vhelium3(j)
            deuterium(j)=vdeuterium(j)
            carbon(j)=vcarbon(j)
            trogen(j)=vtrogen(j)
            oxygen(j)=voxygen(j)
          enddo
C
C  End of variable reset.
C
        else
C
C  CONVERGENCE IN HENYEY.  ADVANCE THE COMPOSITION
C
          call compchange  
          call addmass
          if(konv.lt.2) call addsub
          call printm(imod)
C
C  OPTIMIZE THE TIME STEP, TAKING INTO ACCOUNT THE PREVIOUS CONVERGENCE
C  BEHAVIOR AND THE NET CHANGES.
C
          if(CHANGE.lt.CHGMIN .and. ITER.lt.9 .and. konv.eq.0)
     *                         DTIME=DTIME*1.1
          if(CHANGE.gt.CHGMAX) DTIME=DTIME*0.8
          if(CHANGE.gt.2.*CHGMAX) DTIME=DTIME*0.8
          DTIME=max(DTIME,DTMIN)
          DTIME=min(DTIME,DTMAX)
        endif
        if((nmod.gt.0 .and. imod.ge.nmod) .or. konv.ge.7) goto 99
      enddo
 99   call printm(-1)
      stop
      end
      SUBROUTINE PRINTM(imod)
C+
C  THE SUBROUTINE PRINTM WRITES OUT MODEL INFORMATION
C
C  AUTHOR: H.W. YORKE   30-AUG-02               (JPL / CALTECH)
C-
      include 'parm.h'
      include 'var.h'
      CHARACTER*1 CS(2)
      DATA MODWRT/-1/
      DATA CS/'*',' '/
      save
C
      if(mod(IMOD,NRIT).ne.0 .and. mod(MODEL,JRIT).ne.0) return
      if(IMOD.gt.0) then
        te4 = X(N,3)/7.125e-4/X(N,2)/X(N,2) 
        tee=sqrt(sqrt(te4))
c       Znet=1.d0-hydrogen(1)-helium4(1)-helium3(1)-deuterium(1)        &
c    &      -carbon(1)-trogen(1)-oxygen(1)
        if(iter.lt.abs(itmax/2)) then
          write (6,200)  MODEL,TIME,DTIME,CHANGE,                       &
     &                 X(1,1),X(1,4),X(N,2),X(N,3),tee,Zmass,Zflux 
c    &        hydrogen(1),helium4(1),helium3(1),deuterium(1),           &
c    &        carbon(1),oxygen(1),Znet,                                 &
200       format(' MODEL:',i6,'  TIME:',1p,D11.4,' DTIME:',D11.4,       &
     &     ' NET CHANGE: ',D9.2,                                        &
     &     /,' Pc,Tc,R,L,TE,M:',6D10.3,D11.3                            &
     &     /,1x,83("*"))
c    &     /,' H,He4,He3,D,C,O:',0p,7f8.5,1p,                           &
        else
          write (6,300)  MODEL,TIME,DTIME,CHANGE,                       &
     &                 X(1,1),X(1,4),X(N,2),X(N,3),tee,Zmass,Zflux
300       format(' MODEL:',i6,'  TIME:',1p,D11.4,' DTIME:',D11.4,       &
     &     ' NET CHANGE: ',D9.2,                                        &
     &     /,' XXXXXXXXXXXXXX:',6D10.3,D11.3                            &
     &     /,1x,83("*"))
        endif
      endif 
      if((mod(IMOD,nrit).eq.0 .or. IMOD.lt.0) .and. MODEL.ne.MODWRT)    &
     &  then
        MODWRT=MODEL
        if(deuterium(1).lt.1.d-37) then
          if(helium3(1).lt.1.d-37) then
            write (6,401) MODEL,TIME
          else
            write (6,301) MODEL,TIME
          endif
        else
          write (6,201) MODEL,TIME
        endif
401     format(/,' MODEL:',i6,'  TIME:',1p,D12.4,/                      &
     &         '    J       dM         M         P        R    ',       &
     &         '        L         T         RHO    1-BETA     HYDR  ',  &
     &         '   HE4       C        N        O')
301     format(/,' MODEL:',i6,'  TIME:',1p,D12.4,/                      &
     &         '    J       dM         M         P        R    ',       &
     &         '        L         T         RHO    1-BETA     HYDR  ',  &
     &         '   HE4      HE3       N        O')
201     format(/,' MODEL:',i6,'  TIME:',1p,D12.4,/                      &
     &         '    J       dM         M         P        R    ',       &
     &         '        L         T         RHO    1-BETA     HYDR  ',  &
     &         '   HE4       D        N        O')
        do j=1,N
          call invstate(91,x(j,1),x(j,4),hydrogen(j),helium4(j),Crad,   &
     &         RHO,cP,alpha,beta,delta)
          if(deuterium(1).lt.1.d-37) then
            if(helium3(1).lt.1.d-37) then
            write (6,202) J,CS(1+ICONVECT(j)),dM(j),zM(j),              &
     &                    (X(j,k),k=1,NG),rho,beta,                     &
     &                    hydrogen(j),helium4(j),carbon(j),             &
     &                    trogen(j),oxygen(j)
            else
            write (6,202) J,CS(1+ICONVECT(j)),dM(J),zM(J),              &
     &                    (X(j,k),k=1,NG),rho,beta,                     &
     &                    hydrogen(j),helium4(j),helium3(j),            &
     &                    trogen(j),oxygen(j)
            endif
          else
            write (6,202) J,CS(1+ICONVECT(j)),dM(J),zM(J),              &
     &                    (X(j,k),k=1,NG),rho,beta,                     &
     &                    hydrogen(j),helium4(j),deuterium(j),          &
     &                    trogen(j),oxygen(j)
          endif
C202      format(I5,1x,a1,0p,F9.6,1p,12D9.2)
202       format(I5,1x,a1,1p,D9.2,D13.6,D10.3,D9.2,D13.5,2D10.3,6D9.2)
        enddo
        if(NATM.gt.0 .and. IMOD.ne.0) then
          call atmos(Tatm,RHOatm,Ratm,Patm,X(N,2),X(N,3),NATM)
        else
          call atmos(Tatm,RHOatm,Ratm,Patm,X(N,2),X(N,3),0)
        endif
        if(MODEL.ne.0) then
          write (JUNIT) MODEL,N,TIME,dTIME,taujk,                       &
     &     ((x(j,k),k=1,4),dM(j),                                       &
     &     hydrogen(j),helium3(j),helium4(j),deuterium(j),              &
     &     carbon(j),trogen(j),oxygen(j),                               &
     &     j=1,N)
c         backspace JUNIT
c         read (JUNIT)
          NREC=NREC+1
          write (6,204) MODEL,TIME,NREC,JUNIT
204       format(' MODEL:',i6,'  TIME:',1p,D12.4,' STORED AS RECORD',   &
     &     i5,' ON UNIT',i3)
        endif
      endif
      return
      end
      SUBROUTINE START
C+
C  THE SUBROUTINE START PREPARES FOR THE FIRST RUNNUNG OF HYDRO.
C
C  AUTHOR: H.W. YORKE   30-JAN-02               (JPL / CALTECH)
C-
      include 'parm.h'
      include 'var.h'
      character*50 SFILE,OFILE
      character*1  why(3)
      dimension ZH(14)
      dimension dens(MJ),temp(MJ), xl(MJ), r(MJ) ,pres(MJ)
      equivalence (X(1,1),pres) , (X(1,4),temp)
      equivalence (X(1,2),r) , (X(1,3),xl)
      COMMON/ABUND/XBA(14),H1(14),AH(14)
      DATA ZH     /1.0081451,4.003874,12.0038156,23.,24.32,             &
     & 26.97,28.06,32.07,39.102,40.08,55.85,14.0075257,16.,19.99/
      DATA RSUN/6.96d10/
      save
c
      do i=1,14
        ah(i)=zh(i)
      enddo
c
      call comrd
      READ (5,200) IUNIT,SFILE
c     WRITE(6,200) IUNIT,SFILE
      call comrd
      READ (5,200) JUNIT,OFILE
c     WRITE(6,200) JUNIT,OFILE
c
      if(IUNIT.eq.1) then
        open(unit=IUNIT,file=SFILE, form='formatted' , status='old' ) 
      else
        open(unit=IUNIT,file=SFILE, form='unformatted' , status='old' )
      endif
      if(IUNIT.ne.JUNIT)                                                &
     &  open(unit=JUNIT,file=OFILE, form='unformatted' ,                &
     &  status='unknown' )
C
 200  FORMAT(20X,i3,1X,a)
 201  FORMAT(5(6X,I6))
 202  FORMAT(1P,4(6X,D9.2))
 203  FORMAT(4(6X,D9.2))
C
      call comrd
      READ (5,201) NREC,NMOD,NRIT,ITMIN,ITMAX
c     WRITE(6,201) NREC,NMOD,NRIT,ITMIN,ITMAX
      call comrd
      READ (5,201) JADD,JSUB,NATM
c     WRITE(6,201) JADD,JSUB,NATM
      call comrd
      READ (5,203) Atmx,Atmn,dTAUmx,RLH
c     WRITE(6,202) Atmx,Atmn,dTAUmx,RLH
      taufactor=dTAUmx
      call comrd
      READ (5,203) dLmx,dLmn,dXmx,dXmn
c     WRITE(6,202) dLmx,dLmn,dXmx,dXmn
      call comrd
      READ (5,203) dPmx,dPmn,Crad,Cwrk
c     WRITE(6,202) dPmx,dPmn,Crad,Cwrk
      call comrd
      READ (5,203) dZmax,dZmin,dZdt
      dZmax=dZmax*1.00001
c     WRITE(6,202) dZmax,dZmin,dZdt
      call comrd
      READ (5,203) EPS
c     WRITE(6,202) EPS
      call comrd
      READ (5,203) SMIN
c     WRITE(6,202) SMIN
      call comrd
      READ (5,203) SMAX
c     WRITE(6,202) SMAX
C
      call comrd
      READ (5,203) DTIME,FACTIM,DTMIN,DTMAX
c     WRITE(6,202) DTIME,FACTIM,DTMIN,DTMAX
      call comrd
      READ (5,203) CHGMIN,CHGMAX
c     WRITE(6,202) CHGMIN,CHGMAX
      call comrd
      READ (5,*)   (H1(J), J = 1,4) 
      call comrd
      READ (5,*)   (H1(J), J = 5,8) 
      call comrd
      READ (5,*)   (H1(J), J = 9,12) 
      call comrd
      READ (5,*)   (H1(J), J = 13,14) 
      call comrd
      NG=MH
      JRIT = max(1,NRIT/10)
C
      XH2=H1(2)
      H1(2)=abs(H1(2))
      SUM = 1.
      DO J = 1,14
        SUM = SUM - H1(J)
      END DO
      XX=H1(1)
      YY=SUM
C H1(2) is temporarily used as the primordial deuterium abundance and
C then replaced by the initial helium abundance.
      H1(2)=YY
      DO J = 1,14
        XBA(J) = H1(J) / AH(J)
      END DO
      ZZ= 1. - XX - YY 
      IF(IUNIT.EQ.1) THEN
        NREC=0
        MODEL=0
C SETUP STARTING MODEL
        do j=1,MJ
          read(IUNIT,*,end=34) zM(j),r(j),temp(j),xl(j),dens(j)
          if(r(j).le.0.d0) goto 35
        enddo
        J=MJ+1 
        GOTO 35
  34    backspace(IUNIT)
  35    N=min(J-1,MJ)
        close(IUNIT)
        Zmass=zM(N+1)
C**************
        amu = 2.*XX + 0.75*YY + 0.57*ZZ
        Rglog=log10(83145100.d0*amu)
        do i = 1,N
          r(i) = RSUN*r(i) 
          tlog = log10(temp(i))
          hydrogen(i) = XX
          helium4(i) = YY
          helium3(i) = 1.d-5 
          deuterium(i)=abs(XH2)
          carbon(i)=H1(3)
          trogen(i)=H1(12)
          oxygen(i)=H1(13)
          denlog = log10(dens(i))
          xl(i) = xl(i)*3.85d33 
          Plog = Rglog + denlog + Tlog
          pres(i)=10.d0**Plog
          zM(i)=zM(i)*zMass
        end do 
        dM (1)=zM(1)
        do i=2,N
          dM (i)=zM(i)-zM(i-1)
        enddo
        TIME=0.
        taujk=1.D11
        dTAUmx=taujk*taufactor
        call printm(MODEL)
      ELSE
        if(NREC.lt.0) NREC=99999
        nn=max(1,NREC)
        do imod=1,nn
          read(IUNIT,end=301) MODEL,N,TIME,dTIM,taujk,                  &
     &     ((x(j,k),k=1,4),dM(j),                                       &
     &     hydrogen(j),helium3(j),helium4(j),deuterium(j),              &
     &     carbon(j),trogen(j),oxygen(j),                               &
     &     j=1,N)
        enddo
        IMOD=NREC+1
        goto 302
 301    backspace(IUNIT)
 302    NREC=IMOD-1
        zM(1)=dM(1)
        do i=2,N
          zM(i)=zM(i-1)+dM(i)
        enddo
        if(NREC.eq.0) then
          TIME=0.d0
          MODEL=0
        endif
        write(6,209)  NREC,MODEL,N,TIME,dTIM
 209    format(' RECORD',i4,' MODEL:',i6,' N:',i4,' TIME:',             &
     &           1p,D11.4,' DTIME:',D11.4)
        if(JUNIT.ne.IUNIT) NREC=0
        if(DTIME.lt.0.0d0) DTIME=dTIM
        write(6,*) '... continuing calculations with MODEL=',MODEL
        Zmass=zM(N)
        write(*,*) 'Zmass=', Zmass
        if(Jadd.gt.0) then
          why(1)='I'
          call add(Jadd,why)
        endif
        if(Jsub.gt.0) call sub(Jsub)
        if(FACTIM.lt.0.d0) TIME=0.
      ENDIF
      RETURN
      END
      subroutine comrd
      character*1 line(1)
      character*80 lin
      equivalence (lin,line)
      save
c
1     read (5,203,end=999) lin
      if(line(1).eq.' ') then
      write(6,203)         lin
        backspace 5
        return
      else
        if(line(1).ne.'c' .and. line(1).ne.'C')                         &
     &  write(6,203)  lin
c202    format(1x,a)
      endif
      goto 1
 999  write(6,203) 'COMRD: END OF INPUT DATA'
 203  format(a)
      return
      end
