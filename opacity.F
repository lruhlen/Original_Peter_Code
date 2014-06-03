      subroutine opacity(j,tlog,rholog,bkap)
      include 'parm.h' 
      include 'var.h' 
      parameter(NTVAL=60,NRVAL=14)
      DIMENSION TVAL(NTVAL),RVAL(NRVAL)
      DATA TVAL/
     +2.000,2.097,2.176,2.243,2.301,2.352,2.398,2.439,
     +2.477,2.512,2.544,2.574,2.602,2.628,2.653,2.677,
     +2.699,2.740,2.778,2.813,2.845,2.875,2.903,2.929,
     +2.954,2.978,3.000,3.021,3.041,3.061,3.079,3.097,
     +3.114,3.130,3.146,3.161,3.176,3.204,3.230,3.255,
     +3.279,3.301,3.350,3.400,3.450,3.500,3.550,3.600,
     +3.650,3.700,3.800,3.900,4.000,4.079,4.176,4.301,
     +4.477,4.699,4.845,5.000/
      DATA RVAL/-12.,-11.,-10.,-9.,-8.,-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.
     1/
      data ifirst/0/
      save

c     if(j.lt.100) write(6,*) 'opacity: j,tlog,rhol=',j,tlog,rhol
      if(ifirst.eq.0) then
        ifirst=1

C....read in the low temperature opacity 
        open(unit=52, file='opac.cool', status='old')
        do itemp=1,60
          READ(52,46) (ZKAP(idense,itemp), idense=1,7)
c         WRITE(6,46) (ZKAP(idense,itemp), idense=1,7)
          READ(52,46) (ZKAP(idense,itemp),idense=8,14)
c         WRITE(6,46) (ZKAP(idense,itemp),idense=8,14)
46        FORMAT(7F7.3)
        end do
        close(52)
        open(unit=50,file='hyd.cond',status='old')
        open(unit=51,file='hel.cond',status='old')


c.....read in hydrogen conductive opacities
        read(50,*) (rhyr(i),i=1,23)
        read(50,*) (thyr(i),i=1,28)
        do i=1,28
          read(50,*) (auxchyop(k,i),k=1,23)
        end do
c.....read in helium conductive opacities
        read(51,*) (rhel(i),i=1,25)
        read(51,*) (thel(i),i=1,32)
        do i=1,28
          read(51,*) (auxcheop(k,i),k=1,25)
        end do

        close(50)
        close(51) 
      endif

      YY=helium4(j)+helium3(j)

c.....Start here with logT,logRho, Rference logTs logRhos, Table,
c..    log Opacity is Bkap.
      RHOL=rholog
      IF(TLOG .GE. TVAL(NTVAL)) GO TO 206
      call bracket(TVAL,TLOG,NTVAL,mode,it1,it,iterx)
      RHOL=max( RVAL(1) , RHOL )
      RHOL=min( RVAL(NRVAL)-1.E-9 , RHOL )

      WS1 = TVAL(IT1) - TLOG                                                1900
      WS = TVAL(IT1) - TVAL(IT)                                             2000
      WS1 = WS1 / WS                                                        2100
      WS = 1. - WS1                                                         2200

      call bracket(RVAL,RHOL,NRVAL,mode,iw1,iw,iterx)
      WS2 = RHOL - RVAL(IW1)                                                2800
      WS3 = RVAL(IW) - RVAL(IW1)                                            2900
      WS2 = WS2 / WS3                                                       3000
      WS3 = 1. - WS2                                                        3100
      Z00 = ZKAP(IW,IT)                                                     3200
      Z10 = ZKAP(IW1,IT)                                                    3300
      Z11 = ZKAP(IW1,IT1)                                                   3400
      Z01 = ZKAP(IW,IT1)                                                    3500
      IF(Z00.EQ.0..OR.Z01.EQ.0..OR.Z11.EQ.0..OR.Z10.EQ.0.) GO TO 200        3600
      WS4 = WS2*(WS1*Z00 + WS*Z01)                                          3700
      WS5 = WS3*(WS1*Z10 + WS*Z11)                                          3800
      BKAP= (WS4+WS5)  

C  Slowly transition into OPAL opacities
      if(Tlog.gt.TVAL(NTVAL)-.1) then
        W1=(TVAL(NTVAL)-Tlog)*10.
        W2=1.d0-W1
        call opaltab(Tlog,RHOlog,hydrogen(j),YY,bkap2,ierror)
        if(ierror.eq.0)  BKAP=W1*BKAP+W2*bkap2
      endif
      GO TO 205

 206  continue
      call opaltab(Tlog,RHOlog,hydrogen(j),YY,bkap,ierror)
c     if(ierror.ne.0) Bkap= 4.7

c.....add in conductive opacity:
      if(rhol.gt.0. .and. tlog.gt.5.0)then
        if(hydrogen(j) .le. 1.e-5) then
          hyckap=1.
        else
          call chyopacity(tlog,rhol,hyckap)
        end if
        call cheopacity(tlog,rhol,heckap)
        ckap=hydrogen(j)*hyckap+(1.-hydrogen(j))*heckap
        ckap=10**ckap
        realbkap = 10.**bkap
        realbkap=1./(1./realbkap+1./Ckap)
        bkap=log10(realbkap)
      end if

205   continue                 
      RETURN                   
200   write(6,201) TLOG,RHOL    
201   FORMAT('  OUTSIDE TABLE, T = ', 1PE12.3, '  RHO= ',E12.3)
      STOP 'opacity'
      END

           Subroutine opacityhy(T,rho,realk)
      include 'parm.h' 
           COMMON/AUXhyOPACITY/RHOREF(18),TREF(15),rkapparef(18,15)
c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearly interpolating between the proper points
c..........of a pure hydrogen opacity table.

           data its,irhs /15,18/
      save

           if ((T.le.Tref(1)).or.(T.ge.Tref(its)).or.(Rho.le.Rhoref(1))
     +     .or.(RhO.ge.Rhoref(irhs))) then
                  write(6,*) 'Temperature: ',T,' Density: ',Rho
                  write(6,*) 'Out of Bounds.'
                 realk=1.
C                goto 2104
                 stop 'opacityhy'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature.

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

c.........We have thus determined k.
2104   continue
c       write(*,*)  ihot,icool,jdense,jrare,krare,kdense,k


c       write(*,*) 'rkapp(jr,ic),rk(jr,ih),rk(jd,ic),rk(jd,ih)'
c       write(*,*)  rkapparef(jrare,icool),rkapparef(jrare,ihot)
c       write(*,*)  rkapparef(jdense,icool),rkapparef(jdense,ihot)
c       write(*,*) 'tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)'

c       write(*,*) tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)


        return
        end
 

           Subroutine opacityhe(T,rho,realk)

c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearly interpolating between the proper points
c..........of a pure helium opacity table.
      include 'parm.h' 
           COMMON/AUXheOPACITY/RHOREF(18),TREF(15),rkapparef(18,15)
           data its,irhs /15,18/
      save

           if ((T.le.Tref(1)).or.(T.ge.Tref(its)).or.(Rho.le.Rhoref(1))
     +     .or.(RhO.ge.Rhoref(irhs))) then
                  write(6,*) 'Temperature: ',T,' Density: ',Rho
                  write(6,*) 'Out of Bounds.'
                 realk=1.
C                goto 2104
                 stop 'opacityhe'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature.

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))



c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

c.........We have thus determined k.
2104   continue
c       write(*,*)  ihot,icool,jdense,jrare,krare,kdense,k
c       write(*,*) 'rkapp(jr,ic),rk(jr,ih),rk(jd,ic),rk(jd,ih)'
c       write(*,*)  rkapparef(jrare,icool),rkapparef(jrare,ihot)
c       write(*,*)  rkapparef(jdense,icool),rkapparef(jdense,ihot)
c       write(*,*) 'tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)'

c       write(*,*) tref(ihot),tref(icool),rhoref(jdense),rhoref(jrare)


        return
        end

           Subroutine chyopacity(T,rho,realk)

c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearlly interpolating between the proper points
c..........of the Hubbard-Lampe (gL-padded) hydrogen conductive opacity table

      include 'parm.h' 
           common /ahycopacity/ Rhoref(23),Tref(28),rkapparef(23,28)
           data its,irhs /28,23/
      save

           condmin=rkapparef(1,its)
c          write(*,*) 'chyopacity: Temperature: ',T,' Density: ',Rho
           if (T.le.Tref(1).or.Rho.le.Rhoref(1)) then
             write(*,*) 'Temperature: ',T,' Density: ',Rho
             write(*,*) 'Out of Bounds in conductive hydrogen O-table.'
             stop 'chyopacity'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

          realk=min(realk,condmin)

c.........We have thus determined k.

          return
          end

           Subroutine cheopacity(T,rho,realk)

c..........This subroutine computes the opacity for a specified temperature
c..........and pressure by linearlly interpolating between the proper points
c..........of the Hubbard-Lampe (gL-padded) helium conductive opacity table

      include 'parm.h' 
           common /ahecopacity/ Rhoref(25),Tref(32),rkapparef(25,32)
           data its,irhs /32,25/
      save
c
             condmin=rkapparef(1,its)
           if(T.le.Tref(1).or.Rho.le.Rhoref(1)) then
             write(*,*) 'Temperature: ',T,' Density: ',Rho
             write(*,*) 'Out of Bounds in conductive helium O-table.'
             stop 'cheopacity'
           endif

c..........Determine which of the reference points bracket the temperature.

           call bracket(Tref,T,its,mode,icool,ihot,iterx)
           call bracket(Rhoref,Rho,irhs,mode,jrare,jdense,iterx)

c..........Determine Krare and Kdense by a linear interpolation in temperature

           rare = rkapparef(jrare,icool)
     +             +((rkapparef(jrare,ihot)-rkapparef(jrare,icool))/
     +             (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

           dense = rkapparef (jdense,icool)
     +              +((rkapparef(jdense,ihot)-rkapparef(jdense,icool))/
     +              (Tref(ihot)-Tref(icool)))*(T-Tref(icool))

c.........Determine kappa by a linear interpolation in the density

          realk=rare+((dense-rare)/(rhoref(jdense)-rhoref(jrare)))
     +              *(rho - rhoref(jrare))

          realk=max(realk,condmin)

c.........We have thus determined k.

          return
          end
C+
C  NAME:    bracket (Version 1.0)
C  AUTHOR:  H.W. Yorke (JPL)
C  DATE:    24-Feb-05     V1.0
C  UPDATES:
C
C  The subroutine bracket finds the upper and lower indices of the
C  members of a monotonic increasing table which bracket a given number.
C
C  *** NOTE:    bracket does not check whether the table is monotonically
C               increasing.
C
C  USAGE: call bracket (Xarray , Xval, Narray , mode , I1 , I2 , ITER )
C
C  WHERE    Xarray   is the given monotonically increasing table (INPUT)
C           Xval     is the given number to be bracketed         (INPUT)
C           Narray   is the length of the table                  (INPUT)
C           mode     = -1, 0 , 1 depending on whether Xval is     OUTPUT
C                    less than the smallest value in the table
C                    (mode=-1), greater than the maximum value
C                    (mode=1) or lies within the table (mode=0)
C           I1       is the lower index of bracketing             OUTPUT
C           I2       is the upper index of bracketing             OUTPUT
C           ITER     is the number of iterations                  OUTPUT
C
C  PROGRAMS USED: none
C-
      subroutine bracket(Xarray,Xval,Narray,mode,I1,I2,iter)
      include 'parm.h' 
      dimension Xarray(Narray)
      save
c
      i=0
      mode=0
      if(Xval.lt.Xarray(1)) then
        I1=1
        I2=2
        mode=-1
        return
      endif
      if(Xval.gt.Xarray(Narray)) then
        I1=Narray-1
        I2=Narray
        mode=1
        return
      endif
      I1=1
      I2=Narray
      do i=1,Narray
        ID=(I1+I2)/2
        if(Xarray(ID).lt.Xval) then
          I1=ID
        else
          I2=ID
        endif
        if(I2-I1.eq.1) goto 90
      enddo
      stop 'bracket: Unable to find position in table'
 90   iter=i
      continue
      return
      end

      subroutine opaltab(Tlog,RHOlog,XH,Y,opac,ierror)
      include 'parm.h' 
      include 'var.h' 
      parameter (MX=14,MR=19,MT=70)
      character*5 TABLE
      dimension OPACTB(MX,MR,MT),XX(MX),RR(MR),TT(MT)
      COMMON/ABUND/XBA(14),H1(14),AH(14)
      data ifirst/0/
      data Z0,Z1/0.,1./
      save

      ierror=0
      if(ifirst.eq.0) then
        ifirst=1
        open(51,file='GN93hz',form='formatted',status='old',err=99)
        do i=1,239
          read(51,*)
        enddo
        do ix=1,MX
          read(51,*)
 201      format(a5,31x,f6.4,12x,f6.4,2(5x,f6.4))
          read(51,201,err=98) TABLE,XX(ix),ZZ,XC,XO
c         write(6,201)        TABLE,XX(ix),ZZ,XC,XO
          if(TABLE.ne.'TABLE') goto 98
C
C  XX(ix) is normally used to store the hydrogen mass content used to
C  generate the table.  However, for X=0 (no hydrogen) then -XX(ix) is
C  used to store the amount of additional oxygen and carbon.
C
          if(XC+XO.gt.0.0001) then
            XX(ix)=-(XC+XO)
          else
            ZZtab = ZZ
          endif
          read(51,*)
          read(51,*)
          read(51,*)
 202      format(4x,f6.1,18f7.1)
          read(51,202,err=98) RR
c         write(6,202) RR
          read(51,*)
          do it=1,MT
 203        format(f4.2,19f7.3)
            read(51,203,err=98) TT(it),(OPACTB(ix,ir,it),ir=1,MR)
c           write(6,203) TT(it),(OPACTB(ix,ir,it),ir=1,MR)
          enddo
c         write(6,203) TT
        enddo
        close(51)
        dr=(RR(MR)-RR(1))/float(MR-1)
        IXmax=1
        XXmax=XX(1)
        do ix=1,MX
          if(XX(ix).gt.XXmax) then
            IXmax=ix
            XXmax=XX(ix)
          endif
        enddo
        Znorm = Z1-H1(1)-H1(2)
      endif
 
      if(XH.lt.XX(1) .or. XH.gt.XXmax) ierror=1
      dZZ=((Z1-XH-Y)-Znorm) * (Z1-ZZtab)/(Z1-Znorm)
      if(dZZ.le.0.0014d0) then
        IX1=IXmax-1
        IX2=IXmax
        if(XH.lt.XX(IX1)) then
          IX1=1
          IX2=2
          if(XH.gt.XX(IX2)) then
            IUP=IXmax-1
            IDN=2
            do i=1,2
              IX1=(IUP+IDN)/2
              if(XH.gt.XX(IX1)) then
                IDN=IX1
              else
                IUP=IX1
              endif
           enddo
           IX1=IDN
           IX2=IUP
          endif
        endif
        DX2=(XH-XX(IX1))/(XX(IX2)-XX(IX1))
        DX1=Z1-DX2
      else
        if(XH.gt.1.d-20) then
          ierror=99
          write(6,*) 'You have no business being in this part of the',  &
     &     ' opacity table  ZZtab,dZZ=',ZZtab,dZZ
          write(6,*) 'X,Y,Znorm=',XH,Y,Znorm
          write(6,'(a,1p,14E11.3)') 'Abundances:',H1
          do j=1,N
            if(abs(XH-hydrogen(j)).lt.1.d-10) then
              write(6,'(a,i4,9F8.5)') 'J,lnT,X,Y4,Y3,D,C,N,O=',j,Tlog,  &
     &          XH,helium4(j),helium3(j),deuterium(j),carbon(j),        &
     &          trogen(j),oxygen(j),1.d0-XH-helium4(j)-helium3(j)-      &
     &          deuterium(j)-carbon(j)-trogen(j)-oxygen(j)
            endif
          enddo
          stop 'opaltab'
        endif
        IX1=MX-1
        IX2=MX
        if(dZZ.lt.-XX(IX1)) then
          IX1=1
          IX2=IXmax+1
          if(dZZ.gt.-XX(IX2)) then
            IUP=MX-1
            IDN=IXmax+1
            do i=1,2
              IX1=(IUP+IDN)/2
              if(dZZ.gt.-XX(IX1)) then
                IDN=IX1
              else
                IUP=IX1
              endif
           enddo
           IX1=IDN
           IX2=IUP
          endif
        endif
        DX2=(dZZ+XX(IX1))/(-XX(IX2)+XX(IX1))
        DX1=Z1-DX2
      endif
c     write(6,204) 'X:',XH,XX(ix1),XX(ix2),ix1,ix2,dx1,dx2
 204  format(a,1p,3E12.4,2i10,2E12.4,i3)

      Rlog=RHOlog-3.*(Tlog-6.)
      if(Rlog.lt.RR(1) .or. Rlog.gt.RR(MR)) ierror=ierror+2
      DR2=(Rlog-RR(1))/dr
      IR1=max(DR2,Z0)
      IR2=min(IR1+2,MR)
      IR1=IR2-1
      DR2=DR2-float(IR1-1)
      DR1=Z1-DR2
c     write(6,204) 'R:',Rlog,RR(ir1),RR(ir2),ir1,ir2,dr1,dr2

      if(Tlog.lt.TT(1) .or. Tlog.gt.TT(MT)) ierror=ierror+4
      call bracket(TT,Tlog,MT,mode,IT1,IT2,iterx)
      DT2=(Tlog-TT(IT1))/(TT(IT2)-TT(IT1))
      DT1=Z1-DT2
c     write(6,204) 'T:',Tlog,TT(it1),TT(it2),it1,it2,dt1,dt2,iterx

      opac=dx1*(dr1*(dt1*OPACTB(ix1,ir1,it1)+dt2*OPACTB(ix1,ir1,it2))
     &        + dr2*(dt1*OPACTB(ix1,ir2,it1)+dt2*OPACTB(ix1,ir2,it2)))
     &    +dx2*(dr1*(dt1*OPACTB(ix2,ir1,it1)+dt2*OPACTB(ix2,ir1,it2))
     &        + dr2*(dt1*OPACTB(ix2,ir2,it1)+dt2*OPACTB(ix2,ir2,it2)))

      return
 98   write(6,*) 'Error reading opacity tables'
      stop 'opaltab'
 99   write(6,*) 'Error opening opacity tables'
      stop 'opaltab'
      end
