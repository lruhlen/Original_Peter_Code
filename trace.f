C     PROGRAM TRACE   
C
C     THIS PROGRAM CALCULATES THE SPECTRUM AND INTENSITY MAPS 
C     OF A COOL DUSTY ENVELOPE SURROUNDING A STELLAR BLACK-BODY SOURCE 
C ******
C   NF     NUMBER OF FREQUENCY POINTS
C   NS     NUMBER OF GRID POINTS  IN THE RADIATIVE TRANSFER
C          GRID (S,T,U) ALONG ONE AXIS IN ONE QUADRANT.
C          THE TOTAL NUMBER OF GRID POINTS IN EACH DIRECTION OF THE S-T
C          PLANE DEPENDS ON THE SYMMETRY ASSUMED AND THE VIEWING ANGLE.
C          A MAXIMUM OF 2*NS-1 GRID POINTS IN EACH DIRECTION IS POSSIBLE
C   FACTOR:  LINEAR SIZE OF ONE QUADRANT IN THE S-T PLANE
C ******
C   THETA IS  THE OBSERVATION ANGLE
C ******
C   FQ1,FQ2 ARE THE LOG10 BOUNDARIES OF THE INVESTIGATED FREQUENCY INTERVAL
C ******
C   LTOT   IS THE LUMINOSITY OF THE CENTRAL SOURCE (<-9 TURN OFF)
C   TEFF   IS THE EFFECTIVE TEMPERATURE OF THE CENTRAL SOURCE
C   RSTR   IS THE RADIUS OF THE CENTRAL SOURCE
C ******
C   MZ,MR    ARE THE DIMENSIONS OF THE ENVELOPE MODEL DATA GRID
C ******
C   NP     NUMBER OF LINE OF SIGHT POSITIONS FOR WHICH DETAILED INFO
C          IS WANTED. UP TO 10 ARE POSSIBLE.
C   LS     S,T POSITION OF L-O-S AND THREE FREQUENCIES FOR DETAILED INFO
C   MAPF   IS AN INTEGER ARRAY SPECIFYING THE FREQUENCIES FOR MAP DATA
C ******
C
C   THIS CODE CALCULATES LINE-OF-SIGHT RADIATION TRANSFER IN
C    USER-SPECIFIED GEOMETRIC SYMMETRIES: 1-D, OR 2-D .
C
C  AUTHOR: H.W. YORKE (UNIV. GOETTINGEN)               9-OCT-87
      INCLUDE 'tracedefs.h'
      logical        LSOURC
      CHARACTER*16   IOFILE
      DIMENSION SI(NFQ),FQ(NFQ),LS(6,10),IW(3),TEM(5),TOT(5)
      DIMENSION FN(NFQ) 
      DIMENSION MAPF(NFQ),BEAM(NSS),FLBEAM(NFQ,NSS)
      DATA RADMUL/1.74532925E-2/, TEM/10.,56.,316.,1778.,10000./
      DATA SIGMA/5.66956E-5/
      open (unit=5,file='trace.inp', status='unknown')
      open (unit=6,file='trace.out', status='unknown')
C
C  READ IN INPUT DATA
C  INPUT AND MAP OUTPUT USE STD UNITS
C
C
 200  FORMAT(1X,A)
C
C     Read in name of output files and open
C
C     timeinit=second()
      call comrd
      READ (5,200) IOFILE
      WRITE(6,200) 'OUTPUT FILE FOR TRACE MAP DATA:'
      WRITE(6,200) IOFILE
      OPEN(8,FILE=IOFILE,STATUS='UNKNOWN')
      call comrd
      READ (5,200) IOFILE
      WRITE(6,200) 'OUTPUT FILE FOR TRACE FLUX-LOG:'
      WRITE(6,200) IOFILE
C
C     Read in symmetry parameter and number of grid points    
C
      call comrd
      READ (5,222) ISYM,MZ,MR
      WRITE(6,223) ISYM,MZ,MR
 222  FORMAT(2X,I1,12X,I3,4X,I3)
 223  FORMAT(' (',I1,'-D DATA) MZ=',
     *       I3,' MR=',I3)
      IF(MZ.GT.MZZ .OR. MR.GT.MRR ) THEN
        WRITE(6,509) 
        STOP 'TRACE1'
  509   FORMAT(' TRAP: Model data dimensions too large !!')
      ENDIF
C
C     Read in star parameters
C     STOT is LTOT, luminosity of central star 
      call comrd
      READ (5,220) STOT,TEFF,RSTR
      WRITE(6,221) STOT,TEFF,RSTR
 220  FORMAT(6X,F5.2,8X,F5.2,7X,F6.2)
 221  FORMAT(' LTOT:',F5.2,'   TEFF:',F5.2,'  RSTR:',F6.2)

      STOT=10.**STOT*SOLLUM
      TEFF=10.**TEFF
      RSTR=10.**RSTR
C
C     Read in information for ray tracing procedure:
C     First define the grid
C
      call comrd
      READ (5,201) NS,FACmin,FACTOR
      WRITE(6,202) NS,FACmin,FACTOR
C     FACmin and FACTOR are the minimum and maximum
C     zone sizes of the non-uniform radiative transfer grid 
 201  FORMAT(4X,I3,17X,2E9.2)
 202  FORMAT(' NS=',I3,' min/max FACTORS:',1p,2E9.2)
C
C     Read in theta angle of view
C
      call comrd
      READ (5,203) ITHETA
      WRITE(6,204) ITHETA
      IPHI=0
 203  FORMAT(8X,I3)
 204  FORMAT(' THETA: ',I3)
C
C     The next read defines the (log) frequency grid
      call comrd
      READ (5,205) NF,FQ1,FQ2
      NF=MIN(NF,NFQ)
      WRITE(6,206) NF,FQ1,FQ2
 206  FORMAT(' NF=',I3,' FQ1: ',F8.5,' FQ2: ',F8.5)
 205  FORMAT(4X,I3,6X,F8.5,6X,F8.5)
C
C     Initialize frequencies
C
      DO 5 I=0,NF-1
         FQ(I+1)=FQ1+(FQ2-FQ1)/(NF-1)*I
 5    CONTINUE
C
C     Read in information for getting detailed printouts
C     along particular lines of sight 
      call comrd
      READ (5,211) NP
      WRITE(6,212) NP
 211  FORMAT(7X,I3)
 212  FORMAT('    NP=',I3)
      DO 10 I=1,NP
         call comrd
         READ (5,214) (LS(J,I),J=1,6)
         WRITE(6,213) (LS(J,I),J=1,6)
 10   CONTINUE
 213  FORMAT(' (S,T,G)',3I3,' (F1,F2,F3)',3I3)
 214  FORMAT(8X,3I3,11X,3I3)
C     Give the parameters of the cloud
      call comrd
      READ(5,560) nopt,cloudsize, rhomax, axis
560   format(6x,i2,1p2e12.2,0pf5.2)
      WRITE(6,561) cloudsize, rhomax, axis
561   format(1X, 'CLOUD SIZE ', 1pe12.4, ' RHOMAX ', e12.4,
     & ' AXIS RATIO ', e12.2)
C
C     Give the beam sizes for which spectra are wanted
C
      call comrd
      READ (5,217) Nbeams,(BEAM(i),i=1,Nbeams)
      WRITE(6,216) Nbeams,(BEAM(i),i=1,Nbeams)
 217  FORMAT(8x,I4/(1p5E12.4))
 216  FORMAT(' Nbeams=',I4/(1p5E12.4))
C
C     Give the frequencies for which a detailed map should be constructed
C
      call comrd
      READ (5,219) Nmaps,(MAPF(i),i=1,Nmaps)
      WRITE(6,218) Nmaps,(MAPF(i),i=1,Nmaps)
 219  FORMAT(8x,I4/15I4)
 218  FORMAT(' Nmaps =',I4/15I4)
C
C   READ INPUT MODEL AND PREPARE IT FOR RAY TRACING PROCEDURE
C
      CALL INIT
C
C   READ PREIBISCH ET AL DUST MODEL 
C
      CALL INIDST(FQ,NF)
C
C
C  The background radiation field is assumed to be a blackbody.
C  The dilution factor and characteristic temperature must be 
C  specified in the next card:
C
      call comrd
      READ (5,225) Fbgr,Tbgr
      WRITE(6,226) Fbgr,Tbgr
 225  format(6x,1p,E9.2,6x,E9.2)
 226  format(' Fbgr=',1p,E9.2,' Tbgr=',E9.2)
      CALL PLANKF(Tbgr,FQ,BBOUT,NF)
      do L=1,NF
         BBOUT(L)=BBOUT(L)*Fbgr
      enddo
      call comrd
      READ (5,227) DUzone,RRzone,ETzone,ESzone
      WRITE(6,228) DUzone,RRzone,ETzone,ESzone
 227  format(1p,2(8x,E9.2),0p,2(8x,F7.4))
 228  format(1p,' DUzone=',E9.2,' RRzone=',E9.2,
     *       0p,' ETzone=',F7.4,' ESzone=',F7.4)
C
      SIZmin=SIZE*FACmin
      SIZE=SIZE*FACTOR
      WRITE(6,239) SIZmin,SIZE
 239  format(' SIZE OF RAY TRACING WINDOW:',1p,2E12.4)
      if(Nbeams.le.0) then
        Nbeams=1
        BEAM(1)=SIZE
      endif
C
C   FIX INTERNAL SYMMETRY PARAMETER:
C
C   JSYM=0  THE FLUX DEPENDS ONLY ON THE IMPACT PARAMETER'S DISTANCE
C           FROM THE CENTRAL L-O-S (E.G. AXIAL SYMMETRY POLE-ON)
C   JSYM=1  ONLY ONE QUADRANT IN THE S-T PLANE NEEDS TO BE CALCULATED
C           (E.G. AXIAL SYMMETRY AND EQUATORIAL SYMMETRY EDGE-ON)
C   JSYM=2  ONLY TWO QUADRANTS IN THE S-T PLANE NEED TO BE CALCULATED
C           (E.G. AXIAL SYMMETRY BUT NO EQUATORIAL SYMMETRY)
C   JSYM=3  ALL FOUR QUADRANTS IN THE S-T PLANE NEED TO BE CALCULATED
C
        JSYM=ISYM
        IF(JSYM.LT.3 .AND. ITHETA.EQ.0) JSYM=0
        IF(JSYM.LT.3 .AND. ITHETA.EQ.90) JSYM=1
C
        MS=NS
        MT=NS
        IF(JSYM.EQ.0) MT=1
        IF(JSYM.GE.2) MT=MS+MS-1
        IF(JSYM.EQ.3) MS=MS+MS-1
        IF(MT.GT.NTT .OR. MS.GT.NSS) THEN
          WRITE(6,215) 'TRACE  : TOO MANY GRID POINTS'
  215     FORMAT(1X,A)
          STOP 'TRACE'
        ENDIF
        IF(MOD(NS,2).EQ.0) THEN
          WRITE(6,*) 'TRACE: THE NUMBER OF GRID POINTS MUST BE ODD'
          STOP 'TRACE'
        ENDIF
C
        THETA=FLOAT(ITHETA)*RADMUL
        PHI=FLOAT(IPHI)*RADMUL
C  DECIDE WHETHER CENTRAL SOURCE IS PRESENT
        LSOURC=STOT.gt.1.E-9*SOLLUM .and. RSTR.gt.1.
      DO L=1,NF
C  Zero the integrated flux array FLUX
        FLUX(L)=0.
C  Zero the integrated flux vs. beam size array FLBEAM
        do j=1,NS
          FLBEAM(L,j)=0.
        enddo
      ENDDO
C   THE MATRIX ELEMENTS FOR KRAY ARE FIXED ACCORDING TO
C       (THETA,PHI)
        A(1)= COS(PHI)
        A(4)= SIN(PHI)
        A(8)= SIN(THETA)
        A(9)= COS(THETA)
        A(2)=-A(9)*A(4)
        A(3)= A(8)*A(4)
        A(5)= A(9)*A(1)
        A(6)=-A(8)*A(1)
        A(7)=0.
C  Make sure that 1/cos and 1/sin are finite
        big=1.E11
        small=1./big
        A(10)=sign(big,A(8))
        if(abs(A(8)).gt.small) A(10)=1./A(8)
        A(11)=sign(big,A(9))
        if(abs(A(9)).gt.small) A(11)=1./A(9)
C
C   The origin can be shifted with respect to the model coordinates
C   by specifying ZM
C
        ZM=0.

C  AND HERE ARE THE S-T LOOPS FOR THE INTENSITY

        SGRID(1)=SIZmin
        SGRID(2)=SIZE
        TGRID(1)=SIZmin
        TGRID(2)=SIZE
      DO 51 I=1,MS
      DO 50 J=1,MT
        CALL FIXGRD(I,J,AREA,JSYM)
        SS=SGRID(I)
        TT=TGRID(J)
C  IW IS A CONTROL VARIABLE FOR PRINT OUTPUT IN KRAY
        IW(1)=0
        DO 37 L=1,NP
        IF(LS(1,L).NE.I .OR. LS(2,L).NE.J) GOTO 37
          IW(1)=LS(4,L)
          IW(2)=LS(5,L)
          IW(3)=LS(6,L)
  37    CONTINUE

C  THE STARTING INTENSITY SI MUST BE SPECIFIED
C  YOU MIGHT WANT TO START OFF WITH 3 DEGREE BLACKBODY RADIATION?
        DO  L=1,NF
          SI(L)=0.
        ENDDO
C
C ***  IF LZ = 1 THEN  Z<0  SHOULD CONTRIBUTE TO THE EMISSION
        LZ = 1
C
C  TEST WHETHER THE L-O-S PASSES THROUGH THE CENTRAL SOURCE. IF IT
C  DOES, THEN NORMALIZE THE INTENSITY TO GIVE THE CORRECT LUMINOSITY
C  FOR THE OPTICALLY THIN CASE AND USE IT FOR THE STARTING INTENSITY.
C  START THE L-O-S INTEGRATION AT Z=0.
C
        IF(ABS(SS)+ABS(TT).lt.0.1*SIZmin .and. LSOURC) THEN
          LZ = 0
          CALL PLANKF(TEFF,FQ,SI,NF)
          DO 39 L=1,NF
  39        SI(L)=SI(L)*PI*RSTR**2/AREA
            write(6,659) ss,tt,rstr,teff,stot
  659       format(' A central source is assumed at S=',1pe12.4,
     &        ' T=',e12.4,' Rstar:',e12.4,' Tstar:',e12.4,' Lstar:',
     &        e12.4)
        ENDIF
C
C  This version adds the contribution of a central disk within
C  0 < r < RRzone. The parameter LZ will be changed to LZ=0, if
C  the LOS hits the disk.
C
        RR=SS**2+(TT*A(11))**2
        IF(RR .lt. RRzone**2.and. TR(2,2,1) .gt. 11.) then
          CALL DISK(I,J,FQ,SI,NF,LZ,JSYM)
        ENDIF
C
C  THE SUBROUTINE KRAY DOES THE L-O-S INTEGRATION
C
        CALL KRAY(FQ,SS,TT,SI,E1,E2,S1,S2,NF,IW,LZ)

      DO 40 L=1,NF
        AA(L,I,J)=SI(L)
        FLUX(L)=FLUX(L)+AA(L,I,J)*AREA
  40  CONTINUE
        do L = 1,NF
        if(FLUX(L) .le. 0.) FLUX(L) = 1.d-99 
        end do 
      do ib=1,Nbeams
        if(SS**2+TT**2 .le. BEAM(ib)**2) then
          do L=1,NF
            FLBEAM(L,ib)=FLBEAM(L,ib) + AA(L,I,J)*AREA
          enddo
             do L=1,NF
             if(FLBEAM(L,ib) .le.0.) FLBEAM(L,ib)=1.d-99 
             end do 
        endif
      enddo
  50  CONTINUE
  51  CONTINUE
C
C INTEGRATE OVER THE SPECTRUM TO GET THE TOTAL FLUX FOR THIS DIRECTION
C
        TOTFLX=0.
        FQ2=FQ(1)*.5
      DO 58 L=1,NF
        FQ1=FQ2
        FQ2=FQ(L)
        IF(L.NE.NF) FQ2=(FQ(L)+FQ(L+1))*.5
        TOTFLX=TOTFLX+(FQ2-FQ1)*FLUX(L)
  58  CONTINUE
C
C JUST TO TEST THE ACCURACY OF THE FREQUENCY INTEGRATION I INTEGRATE
C OVER SEVERAL BLACKBODIES AND COMPARE THE RESULT (TOT) WITH THE
C ANALYTIC VALUE.
C
      DO 60 J=1,5
        TOT(J)=0.
        CALL PLANKF(TEM(J),FQ,SI,NF)
        FQ2=FQ(1)*.5
      DO 60 L=1,NF
        FQ1=FQ2
        FQ2=FQ(L)
        IF(L.NE.NF) FQ2=(FQ(L)+FQ(L+1))*.5
        TOT(J)=TOT(J)+(FQ2-FQ1)*SI(L)
  60  CONTINUE
      DO 62 J=1,5
 62     TOT(J)= 1. - SIGMA/PI*TEM(J)**4 / (TOT(J)+1.E-37)
      WRITE(6,888) (TEM(J),TOT(J),J=1,5)
 888  FORMAT(' FREQUENCY INTEGRATION ACCURACY AT VARIOUS TEMPERATURES',
     *    /,1P,5('  T=',E9.2,':',G9.2))
C
C OUTPUT THE RESULT ONTO UNIT 8
C
        write(6,*) ' TOTFLX=',TOTFLX
        L=MAPF(1)
        WRITE(8,*) MS, MT, L,NMAPS 
C       WRITE(8,88) (SGRID(I),I=1,MS),(TGRID(J),J=1,MT)
88      format(1x,1p8e11.3)

      do m=1,Nmaps
         L=MAPF(m)
         if(L.gt.0) then
          do i=1,ms
          do j=1,mt
          if(AA(L,I,J) .le. 0.d0)AA(L,I,J)=1.d-49 
          BB(I,J) = log10(AA(L,I,J))
          end do
          end do
       WRITE(8,88)  (SGRID(I),I=1,MS),(TGRID(J),J=1,MT),
     &(( BB(I,J),I=1,MS),J=1,MT)
       end if 
      enddo
      CLOSE(8)
      do ib=1,Nbeams
        eyeone  = 0.
        eyezero = 1.E-37
        fl2=-37.
        if(flbeam(1,ib).gt.1.E-37) fl2=log10(flbeam(1,ib))-2.
        ql2=log10(fq(1))-1.
        fq2=.1*fq(1)
        do L = 1,nf
        Fn(L)  = fq(L) *cvel   
        xnulnu=fn(L)*4.*3.14159*flux(L)/3.E10
          fq1=fq2
          fq2=fq(L)
          ql1=ql2
          ql2=log10(fq2)
          fl1=fl2
          fl2=-37.
          if(flbeam(L,ib).gt.1.E-37) fl2=log10(flbeam(L,ib))
          fml=.5*(fl2+fl1)
          if(fml .gt. -37.) then
            flav = 10.**fml
            fav  = 10.**(fml + 0.5*(ql2+ql1))
            eyeone  = eyeone  + fav *(fq2-fq1)
            eyezero = eyezero + flav*(fq2-fq1)
         endif 
        end do
851     format(1pe11.4,1x,1Pe11.4,1x,1pE11.4)
C        if(ib.eq.2) then
C        write(6,241) eyeone, eyezero
241     format(1x, 'II ', 1P2E10.4)
C       stop 
C        end if 
        tbol = .375*eyeone/eyezero
        totl = log10(eyezero) + log10(4.*PI/SOLLUM)

        write(6,511) NF,MS,MT,IB,ITHETA,IPHI,JSYM,IMOD,ITIME,ISEQ,
     $         TMASS,TMC,ELC,REQ,TEFF,BEAM(ib),tbol,totl
 511    format(10I8,/,1p,8E10.2)
C       WRITE(6,510) (L,Fn(L),FLBEAM(L,IB),L=1,NF)
        WRITE(6,510) (L,log10(Fn(L)),log10(FLBEAM(L,IB)),L=1,NF)
 510    FORMAT('    L     FQ         FLUX',/,1x,30("="),
     *     200(/I5,1p,2E12.4))
      enddo

      OPEN(8,FILE=IOFILE,STATUS='UNKNOWN',FORM='FORMATTED')
      do i=1,100000
         read(8,200,END=322)
      enddo
 322  backspace 8
      totl = log10(TOTFLX) + log10(4.*PI/SOLLUM)
      write(8,511) NF,MS,MT,NS,ITHETA,IPHI,JSYM,IMOD,ITIME,ISEQ,
     $       TMASS,TMC,ELC,REQ,TEFF,BEAM(Nbeams),tbol,totl
      WRITE(8,510) (L,log10(Fn(L)),log10(FLUX(L)),L=1,NF)
      CLOSE(8)


 999  WRITE(6,100)
 100  FORMAT(' NORMAL STOP IN MAIN PROGRAM  TRACE  ')
C     timefin=second()
C     timeel=timefin-timeinit
C     write(6,1000) timeel
1000  format(1x, 'time=', 1Pe12.4) 
      STOP
      END
      SUBROUTINE DISK(II,JJ,FQ,SI,NF,LZ,JSYM)
      INCLUDE 'tracedefs.h'

      dimension FQ(NF),SI(NF),ZV(MO)
      data gg0,gg1/0.,1./
      data ifirst,Imax,TTmax/0,30,3000./
      save
c
C  Initialize the model disk parameters if first call
      if(ifirst.eq.0) then
        ifirst=1
        XXmax=RRzone/SCALZ(MGG) 
        IRPOS=XXmax
        IGpos=MGG
        if(IRPOS.gt.MR-4) then
          write(6,*) 'RRzone=',RRzone,' is outside of innermost grid'
          stop
        endif
        IR0=max(2,IRPOS+1)
        IR1 = IRPOS+2
        WS1 = XXmax - float(IRPOS)
        WS0 = 1. - WS1
        TTzone=WS0*TR(2,IR0,IGpos)+WS1*TR(2,IR1,IGpos)
        SSzone=0.
        DO IZ=2,MZ-3
          SSzone=SSzone + SCALZ(MGG)*
     *       (WS0*RHOMOD(IZ,IR0,IGpos)+WS1*RHOMOD(IZ,IR1,IGpos))
        ENDDO
        SSzone=SSzone*2.
        XXmin=1.E-37
        if(ETzone.lt.0. .and. TTzone.gt.0.)
     *    XXmin=(TTmax/TTzone)**(1./ETzone)
        write(6,201) RRzone*XXmin,RRzone
 201    format(' An inner disk has been set up for radii between',
     *     1pE12.4,' and',E12.4)
        write(6,202) TTzone*XXmin**ETzone,TTzone,
     *     SSzone*XXmin**ESzone,SSzone
 202    format(' Here, the temperature is between:',1pE12.4,
     *      ' and',E12.4,/,
     *      '   The surface density is between:',E12.4,
     *      ' and',E12.4)
      endif
C  End of initialization. Now get ready to integrate over area of
C  grid cell. First, get the boundaries of the cell.
      if(II.eq.1) then
        Smin=0.
      else
        Smin=.5*(SGRID(II-1)+SGRID(II))
      endif
      if(II.eq.MS) then
        Smax=SGRID(MS)
      else
        Smax=.5*(SGRID(II+1)+SGRID(II))
      endif
      if(JJ.eq.1) then
        Tmin=0.
      else
        Tmin=.5*(TGRID(JJ-1)+TGRID(JJ))
      endif
      if(JJ.eq.MT) then
        Tmax=TGRID(MT)
      else
        Tmax=.5*(TGRID(JJ+1)+TGRID(JJ))
      endif
      Jmax=Imax
      if(JSYM.eq.0) Jmax=1

C  Now integrate the disk's contribution to the flux
C  over the entire cell.
      do i=1,Imax
        Sp=Smin+(Smax-Smin)*(float(i)-.5)/float(Imax)
        weight=1./float(Imax**2)
        if(JSYM.eq.0) weight=float(2*i-1)*weight
        do j=1,Jmax
          Tp=Tmin+(Tmax-Tmin)*(float(j)-.5)/float(Jmax)
          XXzone=sqrt(Sp**2 + (A(11)*Tp)**2)/RRzone
          if(XXzone.gt.XXmin .and. XXzone.le.1.) then
            if(LZ.eq.1) LZ=0
            Tdisk=TTzone*XXzone**ETzone
            Sdisk=A(11)*SSzone*XXzone**ESzone
            CALL PLANKF(Tdisk,FQ,S1,NF)
C
            DO ID=1,NDST
C  Dust sublimation:
              ZV(ID)=MAX(gg0,(gg1-Tdisk/TSUB(ID))*20.)
              ZV(ID)=MIN(gg1,ZV(ID))
            ENDDO
C  There should be at least some opacity, even when the
C  most refractory particles are destroyed.
            ZV(3)=MIN(gg1,ZV(3)+1.E-4)

            DO L=1,NF
              DEXT=0.
              DO ID=1,NDST
                DEN=(GM(ID)+GM(1)*(gg1-ZV(1))*ERH(ID))*AB(ID)*Sdisk
C  Contributions of the individual dust types to the extinction:
                DEXT=DEXT+ZV(ID)*FEE(L,ID)*DEN
C  Contributions of the individual dust types to the emission:
              ENDDO
              SI(L)=SI(L)+S1(L)*(1.-exp(-DEXT))*weight
c             write(6,203) i,j,XXzone,Tdisk,Sdisk,L,FQ(L),DEXT
 203          format(2I5,1p,3e10.2,i5,2E10.2)
            ENDDO
          endif
        enddo
      enddo
      RETURN
      END
      SUBROUTINE FIXGRD(II,JJ,AREA,JSYM)
C+
C-
      INCLUDE 'tracedefs.h'

      save
C
C   INTERNAL SYMMETRY PARAMETER:
C
C   JSYM=0  THE FLUX DEPENDS ONLY ON THE IMPACT PARAMETER'S DISTANCE
C           FROM THE CENTRAL L-O-S (E.G. AXIAL SYMMETRY POLE-ON)
C   JSYM=1  ONLY ONE QUADRANT IN THE S-T PLANE NEEDS TO BE CALCULATED
C           (E.G. AXIAL SYMMETRY AND EQUATORIAL SYMMETRY EDGE-ON)
C   JSYM=2  ONLY TWO QUADRANTS IN THE S-T PLANE NEED TO BE CALCULATED
C           (E.G. AXIAL SYMMETRY BUT NO EQUATORIAL SYMMETRY)
C   JSYM=3  ALL FOUR QUADRANTS IN THE S-T PLANE NEED TO BE CALCULATED
C
      DATA IFIRST/0/
      save
C
      IF(IFIRST.EQ.0) THEN
        IFIRST=1
        isym=0
        if(jsym.eq.3) isym=1
        call setgrd(SGRID,MS,isym)
c       write(6,*) 'setgrd done'
        if(jsym.eq.2) isym=1
        call setgrd(TGRID,MT,isym)
c       write(6,*) 'setgrd done'
      ENDIF
      if(II.eq.1) then
        Smin=SGRID(1)
      else
        Smin=.5*(SGRID(II)+SGRID(II-1))
      endif
      if(II.eq.MS) then
        Smax=SGRID(MS)
      else
        Smax=.5*(SGRID(II)+SGRID(II+1))
      endif
      if(JSYM.eq.0) then
        AREA=PI*(Smax**2-Smin**2)
        GOTO 999
      endif
      if(JJ.eq.1) then
        Tmin=TGRID(1)
      else
        Tmin=.5*(TGRID(JJ)+TGRID(JJ-1))
      endif
      if(JJ.eq.MT) then
        Tmax=TGRID(MT)
      else
        Tmax=.5*(TGRID(JJ)+TGRID(JJ+1))
      endif
      AREA=(Smax-Smin)*(Tmax-Tmin)
      if(JSYM.eq.1) AREA=AREA*4.
      if(JSYM.eq.2) AREA=AREA*2.
C
 999  RETURN
      END
      FUNCTION GFAC(DR,NDR,RL,RU)
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      Z0=0.
      IF(NDR.LE.0) THEN
        GFAC=Z0
        RETURN
      END IF
      DR2=DR
      RNR=RU-RL
      IF(RNR.EQ.Z0) THEN
      WRITE(6,181)
  181 FORMAT(1X/,'  WRONG INPUT PARAMETERS FOR GRID SETUP')
      STOP
      END IF
      FR=1.1
      IT=1
   25 CONTINUE
      IF(IT.LT.200) GOTO 26
      WRITE(6,180)
  180 FORMAT(1X/,'  ITERATIONS IN GFAC DO NOT CONVERGE')
      STOP
   26 CONTINUE
      FRP=FR**NDR
      FR1=FR-1.
      DEL=DR2*(FRP-1.)/FR1-RNR
      DDEL=DR2*(NDR*FRP/FR-(FRP-1.)/FR1)/FR1
      DFR=-DEL/DDEL
      IF(ABS(DFR).GT..05D0) DFR=SIGN(.05D0,DFR)
      FR=FR+DFR
      IT=IT+1
      IF(ABS(DFR).GT.1.E-6) GOTO 25
      GFAC=FR
      RETURN
      END
      SUBROUTINE SETGRD(GRID,NS,JSYM)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      dimension grid(NS)

C
      write(6,*) ' Grid setup for  N=',NS,JSYM
      Z0=0.
      if(NS.le.1) then
        GRID(1)=0.
        NS=1
      else
        Rmin=grid(1)
        Rmax=grid(2)
        N=NS-1
        if(JSYM.eq.1) N=NS/2
        FAC=GFAC(Rmin,N,Z0,Rmax)
        DS2=Rmin
        GRID(1)=Z0
        DO IS=2,NS
          GRID(IS)=GRID(IS-1)+DS2
          DS2=DS2*FAC
        ENDDO
        if(JSYM.eq.1) then
          do IS=NS,N+1,-1
            GRID(IS)=GRID(IS-N)
          enddo
          do IS=1,N
            GRID(IS)=-GRID(NS+1-IS)
          enddo
        endif
      endif
      write(6,202) (is,grid(is),is=1,NS)
 202  format(1p,10(i4,e9.2))
      return
      end

      SUBROUTINE INIDST(FQ,NF)
    
      INCLUDE 'tracedefs.h'

      DIMENSION FQ(NF),FL(MF),FE(MF),AL(MF)

      DATA Z10,JFz/10.,0/
C     
      save
C     
C     READ IN TOTAL NUMBER OF DUST COMPONENTS AND NUMBER OF FREQUENCIES
C     IN DATA ARRAY. CHECK IF LESS THAN DIMENSION OF ARRAYS
C     
      call comrd
      READ (5,280) NDST,Nff,UNITL
 280  FORMAT(5X,I3,5X,I3,7X,F4.1)
      NFz=min(Nff,sign(Nff,NFz))
      NFRQ=abs(NFz)
      WRITE(6,480) NDST,NFz,UNITL
 480  FORMAT(' NDST',I3,' NFRQ',I3,' UNITL ',F4.1)
      IF(NDST.GT.MO .OR. NFRQ.GT.MF) THEN
         WRITE(6,489) MO,MF
 489     FORMAT(' *** fatal error in INIDST ***',//,
     &        ' The DIMENSION of several arrays is limited to NDST=',I3,
     &        ' and NFRQ=',i3)
         STOP 'INIDST'
      ENDIF
C     
C     READ IN PARAMETERS OF THE DIFFERENT DUST GRAIN TYPES
      DO 10 ID=1,NDST
         call comrd
         READ (5,202) AB(ID),GM(ID),TSUB(ID),ERH(ID)
         WRITE(6,202) AB(ID),GM(ID),TSUB(ID),ERH(ID)
         GM(ID)=Z10**GM(ID)
         AB(ID)=Z10**(AB(ID)*2.+UNITL)*PI
         TSUB(ID)=Z10**TSUB(ID)
 10   CONTINUE
 202  FORMAT(6(5X,F7.3))
C
C     READ IN FREQUENCIES OF DUST OPACITIES AND CHANGE TO LOG10
      call comrd
      READ (5,482) (FL(L),L=1,NFRQ)
      WRITE(6,482) (FL(L),L=1,NFRQ)
      DO 491 L=1,NFRQ
 491     FL(L)=LOG10(FL(L))
      DO 498 ID=1,NDST
C     
C     NOW READ IN EXTINCTION EFFICIENCY AND ALBEDO FOR THE DUST TYPES
C     
         call comrd
         READ (5,482) (FE(L),L=1,NFRQ)
         WRITE(6,482) (FE(L),L=1,NFRQ)
         call comrd
         READ (5,482) (AL(L),L=1,NFRQ)
         WRITE(6,482) (AL(L),L=1,NFRQ)
C     
C     CALCULATE CONVERSION FACTOR FROM DIMENSIONLESS EXTINCTION
C     COEFFICIENT 'Qext' TO GRAM-OPACITY [cm**2 / g]
C     
C     CALCULATE ABSORPTION GRAM OPACITIES
 482     FORMAT(1X,1P,5E12.5)
         DO 493 L=1,NFRQ
            AL(L)=LOG10(AL(L))
 493        FE(L)=LOG10(FE(L))

C     
C     LOG-LOG INTERPOLATION ONTO FQ-FREQUENCY GRID
         DO 497 L=1,NF
            DO 492 I=2,NFRQ
 492           IF(FQ(L).LT.FL(I)) GOTO 494
            FEE(L,ID)=Z10**FE(NFRQ)
            FEA(L,ID)=Z10**AL(NFRQ)
            GOTO 496
 494        KF2=I
            KF1=KF2-1
            DQ=(FQ(L)-FL(KF1))/(FL(KF2)-FL(KF1))
            FEE(L,ID)=Z10**(FE(KF1)+(FE(KF2)-FE(KF1))*DQ)
            FEA(L,ID)=Z10**(AL(KF1)+(AL(KF2)-AL(KF1))*DQ)
 496        CONTINUE
 497     CONTINUE
 498  CONTINUE
      DO 499 L=1,NF
 499     FQ(L)=Z10**FQ(L)
      WRITE(6,399) 'Dust opacities initialized for ',NF,' frequencies:',
     &        (FQ(L),L=1,NF)
      DO 502 ID=1,NDST
         WRITE(6,399) 'Extinction Cross Sections for component #',
     &        ID,':',(FEE(L,ID),L=1,NF)
         WRITE(6,399) 'Albedo for component #',
     $        ID,':',(FEA(L,ID),L=1,NF)
 502  CONTINUE
 399  FORMAT(/,A,I2,A,1P,100(/,1X,10E12.4))
C   initialize the scattering contribution
C     if(NFz.gt.0) then
C        read(3,ERR=999) NFz,IMODz,ITIMEz,ISEQz,
C    $           (FQz(i),i=1,NFz)
C        if(IMOD.ne.IMODz .or. ISEQ.ne.ISEQz) then
C           write(6,*) 'INIDST: mismatch of scattering model data!'
C           write(6,*) 'IMODz,ITIMEz,ISEQz=',IMODz,ITIMEz,ISEQz
C           write(6,*) 'IMOD,ITIME,ISEQ=',IMOD,ITIME,ISEQ
C           stop
C        endif
         NFz=NF
         do i=1,nfz
         FQz(i) = fq(i)
         end do 
           write(6,*) ' read unit3'  
         do l=1,MFz
C           read(3) (((SMz(i,j,k,l),i=1,MZZ),j=1,MRR),k=1,MGG)
C           do k=1,MGG
C               do i=2,MZZ
C                SMz(i,1,k,l)=SMz(i,2,k,l)
C              enddo
C              do j=2,MRR
C                SMz(1,j,k,l)=SMz(2,j,k,l)
C              enddo
C           enddo
C        enddo

            do i=1,MZZ
            do j=1,MRR
            SMz(i,j,1,l) = 0.
            end do
            end do
         end do 
         last=1
         JFz=0
         do j=1,MFz
            LzMAP(j)=0
            do l=last,NF
               if(abs(1.-FQ(l)/FQz(j)).lt.1.e-5) then
                  JFz=JFz+1
                  LzMAP(JFz)=l
                  last=l+1
                  goto 2
               endif
            enddo
  2         continue
         enddo
         if(JFz.gt.0) then
            write(6,223) JFz,' overlapping frequencies:',
     *                   (i,LzMAP(i),i=1,JFz)
 223        format(1x,i4,a,90(/,8i10))
         else
            write(6,224) ' no overlapping frequencies:',
     *                   (i,FQz(i),i=1,NFz),
     *                   (i,FQ(i),i=1,NF)
 224        format(a,90(/,i5,1pe12.4))
         endif
 999     NFz=JFz
C     endif
      RETURN
      END
C   SET UP ENVELOPE STRUCTURE         
      SUBROUTINE INIT


      INCLUDE 'tracedefs.h'
      PARAMETER (NWR=15,RHOFAC=1.0,SFAC=1.0,TFAC=1.0,DFAC=1.E-0)

      DIMENSION Z(MZZ,1), R(MRR,1), ZHF(MZZ), RHF(MRR)

C

      DATA TTMAX/1.e4/

      NZ=MZ
      NZ1=NZ-1
      NZ2=NZ-2
      NG=1
      NR=MR
      NR1=NR-1
      NR2=NR-2
      ZINT=CLOUDSIZE/DFLOAT(NZ2)
      Z(2,1) = 0.
      DO 314 I=3,NZ
314    Z(I,1) = Z(I-1,1) + ZINT
      Z(1,1) = - Z(3,1)
      R(2,1) = 0.
      DO 315 I=3,NR
315    R(I,1) = R(I-1,1) + ZINT
      R(1,1) = - R(3,1)
      DO 317 I = 1,NR1
317   RHF(I) = 0.5 * (R(I,1) + R(I+1,1))
      RHF(NR) = RHF(NR1)+ R(3,1)
      DO 318 J = 1,NZ1
318   ZHF(J) = 0.5 * (Z(J,1) + Z(J+1,1))
      ZHF(NZ) = ZHF(NZ1)+ Z(3,1)
      if(STOT.lt.1.E-8*SOLLUM) then
         STOT=1.E-15*SOLLUM
         REQ=SOLRAD*.01
         ELC=STOT
         RSTR=REQ
         TEFF=SOLTEF*sqrt(SOLRAD/(REQ+1.e-37)*sqrt(ELC/SOLLUM))
      endif

      WRITE(6,*) 'Parameters of central source:'
      WRITE(6,224) STOT/SOLLUM,TEFF,RSTR/SOLRAD
 224  format(' Luminosity [Lsun] =',F10.2,' Teff [K] =',F10.2,
     *      ' Radius [Rsun] =',F10.2)
       rcon=rhomax*sqrt(rhf(2)**2 + zhf(2)**2)
       xlcon=stot/16./3.14159/5.67d-5
       if(nopt.eq.1)then
       do i=2,nr
       do j=2,nz
       dist=sqrt(rhf(i)**2 + zhf(j)**2)
       t4 = xlcon/dist/dist
       tr(j,i,1) =  sqrt(sqrt(t4))
       rhomod(j,i,1) = rcon/dist
       emod(j,i,1) = 1.25 * rhomod(j,i,1)*boltz/hmas
     & * sqrt(sqrt(t4))
       end do
       end do
       else
       aay=cloudsize
       bby=aay/axis
       ee=sqrt(1. - (bby/aay)**2)
       do i=2,nr
       do j=2,nz
       xtane = zhf(j)/rhf(i)/sqrt(1.-ee*ee)
        eangl=atan(xtane)
       aax=rhf(i)/cos(eangl)
       if(aax.gt.aay)then
       rhomod(j,i,1) = rhomax/1.e4
       else
       rhomod(j,i,1)= rcon/aax
        end if 
       dist=sqrt(rhf(i)**2 + zhf(j)**2)
       t4 = xlcon/dist/dist
       if(t4 .lt. 1.E4) t4 = 1.E4
       tr(j,i,1) = sqrt(sqrt(t4))
       emod(j,i,1) = 1.25 * rhomod(j,i,1)*boltz/hmas
     & * sqrt(sqrt(t4))
       end do
       end do
       end if
C      plot data for model rho and tr
        open(10,file='plotout', status='unknown')
         WRITE(10,*) NZ-1,NR-1
        WRITE(10,41) (ZHF(IZ), IZ=2,NZ), 
     & (RHF(IR), IR=2,NR), ((DLOG10(RHOMOD(IZ,IR,1)), IZ=2,NZ), 
     & IR=2,NR),  ((DLOG10(TR(IZ,IR,1)), IZ=2,NZ ), 
     & IR=2,NR)
41      format(1x,1P6E12.4)  
        close(10)
      DO 10 IG=1,MGG
C clear out the center zones
        EMODc=EMOD(3,2,IG)/RHOMOD(3,2,IG)
        RHOMOD(3,2,IG)=min(RHOMOD(3,2,IG),RHOMOD(4,2,IG)*2.047529)
        RHOMOD(2,2,IG)=min(RHOMOD(2,3,IG),RHOMOD(3,2,IG),RHOMOD(2,2,IG))
        EMOD(2,2,IG)=EMODc*RHOMOD(2,2,IG)
        EMOD(3,2,IG)=EMODc*RHOMOD(3,2,IG)
        do k=2,MZ
          RHOMOD(K,1,IG)=RHOMOD(K,2,IG)
          EMOD(K,1,IG)=EMOD(K,2,IG)
        enddo
        do j=1,MR
          RHOMOD(1,J,IG)=RHOMOD(2,J,IG)
          EMOD(1,J,IG)=EMOD(2,J,IG)
        enddo

C
C In the following the ionization is turned off completely.
C
         DO 13 J=1,MZ
            DO 14 K=1,MR
               XXMOD(J,K,IG)=0.
               EMOD(J,K,IG)=EMOD(J,K,IG)/RHOMOD(J,K,IG)
 14         CONTINUE
 13      CONTINUE

         CALL STATEM(EMOD(1,1,IG),RHOMOD(1,1,IG),XXMOD(1,1,IG),
     $        PPMOD(1,1,IG),YYMOD(1,1,IG),TTMOD(1,1,IG),EEMOD(1,1,IG))
         SCALZ(IG)=Z(3,IG)*SFAC
         DO 15 J=1,MZ
            DO 16 K=1,MR
               TTMOD(J,K,IG)=MIN(TTMOD(J,K,IG),TTMAX)
 16         CONTINUE
 15      CONTINUE

 10   CONTINUE

C
C SIZE is the size of the  Ray Tracing window.
C

      SIZE=R(MR-3,1)*SFAC

 4403 CONTINUE
      RETURN
      END
C
C+++ KRAY
      SUBROUTINE KRAY(FQ,S,T,SI,E1,E2,S1,S2,NF,IW,LZ)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      PARAMETER( NSS=121,NTT=241,NFQ=48)
      COMMON /MATRIX/A(11),ZM,FTR
      DIMENSION S1(NFQ), S2(NFQ), E1(NFQ), E2(NFQ)
C
C  THE SUBROUTINE KRAY CALCULATES RADIATION TRANSFER ALONG THE LINE OF
C  SIGHT (LOS). THE ROTATION MATRIX OF THE LOS IS GIVEN IN THE COMMON
C  BLOCK  MATRIX. THE EXTINCTION COEFFICIENT AND THE SOURCE FUNCTION FOR
C  EACH POINT (X,Y,Z) ALONG THE LOS IS GIVEN IN THE SUBROUTINE SOURCE.
C
C  E1,E2,S1,S2 ARE VECTORS OF LENGTH > NF
C
C  DON'T FORGET TO SPECIFY A STARTING VALUE FOR SI BEFORE CALLING
C  KRAY!
C
      DIMENSION FQ(NF),SI(NF),
     * IW(3),XPL(12)
      DATA IMAX,XPLMIN,ZK/10000,1.024E-34,1.E-37/
      data icount/-1/
C
      save
C
      IF(IW(1).GE.1) WRITE(6,100) S,T,FQ(IW(1)),FQ(IW(2)),FQ(IW(3))
  100 FORMAT(1H1,' LINE OF SIGHT INFORMATION FOR  S=',1PE12.4,', T=',
     *     E12.4,' FQ=',E12.4,2(9X,E12.4)//
     *     4X,'I',5X,'X',9X,'Y',9X,'Z',9X,'U',
     *     7X,'RL',5X,'TL',5X,'XL',1X,3(4X,'DTL',4X,'SNL',4X,'IL+'))
        X0=A(1)*S+A(2)*T
        Y0=A(4)*S+A(5)*T
        Z0=A(7)*S+A(8)*T + ZM
C
C  Initialize the starting point of the LOS
C
      IVORB=LZ
      CALL INITLIN(X0,Y0,Z0,U1,IVORB)
      IF(IVORB.EQ.999) GOTO 1001

        X=X0+A(3)*U1
        Y=Y0+A(6)*U1
        Z=Z0+A(9)*U1
C
C  Get the source function and extinction coefficient at the starting point
C
      CALL SOURCE(FQ,E2,S2,U2,RHX,TEX,
     $       XEX,NF,IENDE)
      if(icount.ge.0) then
        icount=icount+1
        if(icount.gt.400) stop
        write(6,666) RHX,TEX,XEX,NF,IENDE
 666    format('rh,t,x,nf,iende:',1p,3e10.2,2i5)
      endif
      IF(IENDE.EQ.1) GOTO 1001

      XX2=MAX(XEX,ZK)
C
C  Detailed information for the present LOS
C
      IF(IW(1).LE.0) GOTO 10
        XPL(1)=LOG10(MAX(RHX,XPLMIN))
        XPL(2)=LOG10(MAX(TEX,XPLMIN))
        XPL(3)=LOG10(MAX(XEX,XPLMIN))
        DO 9 J=1,3
          L=IW(J)
c         XPL(1+3*J)=LOG10(MAX(E2(L)*DU,XPLMIN))
          XPL(1+3*J)=LOG10(XPLMIN)
          XPL(2+3*J)=LOG10(MAX(S2(L),XPLMIN))
          XPL(3+3*J)=LOG10(MAX(SI(L),XPLMIN))
   9    CONTINUE
        I=0
        WRITE(6,101) I,X,Y,Z,U2,XPL

C  START OF LOS-INTEGRATION ***************

  10    DO 1000 I=1,IMAX
          XX1=XX2
          DO 14 L=1,NF
            E1(L)=E2(L)
            S1(L)=S2(L)
  14      CONTINUE

          DU=U2-U1
          U1=U2

          X=X0+A(3)*U2
          Y=Y0+A(6)*U2
          Z=Z0+A(9)*U2
C
C   THE POINT (X,Y,Z) IS KNOWN. NOW GET THE EXTINCTION COEFFICIENT AND
C   THE SOURCE FUNCTION FROM THE SUBROUTINE  SOURCE.

          CALL SOURCE(FQ,E2,S2,U2,RHX,TEX,
     $         XEX,NF,IENDE)

          XX2=MAX(XEX,XPLMIN)

          DO 40 L=1,NF

C
C  Test if ionization front 
C
             IEFL=1
             IF (ABS(XX2-XX1).LT.0.1) THEN
                IEFL=0
                EMED=(E1(L)+E2(L))/2
                GOTO 7998
             ENDIF
C  If at ionization front, then divide up the integration interval.
C
             DTAU=DU/2.*E1(L)
             IF(DTAU.LE.2.E-3) THEN
                IF(DTAU.GT.0.) THEN
                   ETAU=DTAU*(.5-DTAU/6.)
                   ETAU1=DTAU*(ETAU-1.)
                   SI(L)=SI(L)+(SI(L)-S1(L))*ETAU1
                ENDIF
             ELSE
                IF(DTAU.LT.30.) THEN
                   ETAU=EXP(-DTAU)
                   ETAU1=(1.E0-ETAU)/DTAU
                   SI(L)=S1(L)+(SI(L)-S1(L))*ETAU
                ELSE
                   SI(L)=S1(L)
                ENDIF
             ENDIF
           DTAU=DU/2.*E2(L)
          IF(DTAU.LE.2.E-3) THEN
            IF(DTAU.GT.0.) THEN
              ETAU=DTAU*(.5-DTAU/6.)
              ETAU1=DTAU*(ETAU-1.)
              SI(L)=SI(L)+(SI(L)-S2(L))*ETAU1
            ENDIF
          ELSE
            IF(DTAU.LT.30.) THEN
              ETAU=EXP(-DTAU)
              ETAU1=(1.E0-ETAU)/DTAU
              SI(L)=S2(L)+(SI(L)-S2(L))*ETAU
            ELSE
               SI(L)=S2(L)
            ENDIF
          ENDIF
          GOTO 39

C
C  No, not at ionisation front:
C
 7998       CONTINUE
            DTAU=DU*EMED

          IF(DTAU.LE.2.E-3) THEN
            IF(DTAU.GT.0.) THEN
              ETAU=DTAU*(.5-DTAU/6.)
              ETAU1=DTAU*(ETAU-1.)
              SI(L)=SI(L)+(SI(L)-S1(L))*ETAU1
     *             + (S2(L)-S1(L))*ETAU
            ENDIF
          ELSE
            IF(DTAU.LT.30.) THEN
              ETAU=EXP(-DTAU)
              ETAU1=(1.E0-ETAU)/DTAU
              SI(L)=S2(L)+(SI(L)-S1(L))*ETAU
     *             + (S1(L)-S2(L))*ETAU1
            ELSE
               SI(L)=S2(L)+(S1(L)-S2(L))/DTAU
            ENDIF
          ENDIF
 39       CONTINUE
          S1(L)=DTAU
  40      CONTINUE
          IF(IW(1).LE.0) GOTO 42
            XPL(1)=LOG10(MAX(RHX,XPLMIN))
            XPL(2)=LOG10(MAX(TEX,XPLMIN))
            XPL(3)=LOG10(MAX(XEX,XPLMIN))
          DO 41 J=1,3
            L=IW(J)
            XPL(1+3*J)=LOG10(MAX(S1(L),XPLMIN))
            XPL(2+3*J)=LOG10(MAX(S2(L),XPLMIN))
            XPL(3+3*J)=LOG10(MAX(SI(L),XPLMIN))
  41      CONTINUE
          WRITE(6,101) I,X,Y,Z,U2,XPL,IGG
  101 FORMAT(1X,I4,1P4E10.2,0P12F7.2,2X,I1)
      IF (IEFL.EQ.1) WRITE(6,*) 'E1=',E1(NF),' E2=',E2(NF)
C
C   END OF RADIATION TRANSFER FROM U1 TO U2.
C   Are we at the edge of the model data?

  42    IF(IENDE.EQ.1) GO TO 1001
 1000 CONTINUE
      WRITE(6,999) 'MAXIMUM STEP NUMBER IN KRAY'
  999 FORMAT(1X,A)
      STOP 'KRAY'
 1001 RETURN
      END
C+++ PLANCKF
      SUBROUTINE PLANKF(T,FQ,BF,NF)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      DIMENSION FQ(NF),BF(NF)
C
      save
C
      DO 20 L=1,NF
        XX=1.438833*FQ(L)/T
      IF(XX.GE.0.1) THEN
        IF(XX.GE.50.) THEN
          BF(L)=0.
C         BFT(L)=0.
        ELSE
          EX=EXP(-XX)
          BF(L)=1.191066076E-5*FQ(L)**3*EX/(1.-EX)
C         BFT(L)=BF(L)*XX/T/(1.-EX)
        ENDIF
      ELSE
        EX=1.-XX*(.5-XX*(8.33333333333333333333E-2
     *    -XX*XX*1.38888888888888888889E-3))
        BF(L)=8.278002497E-6*FQ(L)*FQ(L)*T*EX
C       BFT(L)=BF(L)*(XX+EX)/T
      ENDIF
  20  CONTINUE
      RETURN
      END


      subroutine comrd
      character*1 line(1)
      character*78 lin
      equivalence (lin,line)
c
1     read (5,203,end=999) lin
      if(line(1).eq.' ') then
        backspace 5
        return
      else
        write(6,202)  lin
 202    format(1x,a)
      endif
      goto 1
 999  write(6,203) ' END OF INPUT DATA'
 203  format(a)
      stop
      end
C+++  SOURCE
      SUBROUTINE SOURCE(FQ,EXT,SNU,U2,RHOACT,TDST,XACT,NF,IENDE)

      INCLUDE 'tracedefs.h'
C*****
      DIMENSION FQ(NF),BF(NFQ),EXT(NFQ),EXT1(NFQ),SNU(NFQ),
     $     ZV(MO)
      DATA gg0,gg1,ggK,icount/0.,1.,1.E-37,0/
C
      save
C
C     This version calculates source (SNU) and extinction (EXT)
C     function for free-free continuum radiation and dust radiation.
C     The contribution of scattered radiation is included.
C     Input : FQ,X,Y,Z,NF,K
C     Output : EXT,SNU,DLX-Z,RH,TMP,XR
C

      RHOACT= Z00*RHOMOD(IZ0,IR0,IGPOS) + Z01*RHOMOD(IZ0,IR1,IGPOS)
     &      + Z10*RHOMOD(IZ1,IR0,IGPOS) + Z11*RHOMOD(IZ1,IR1,IGPOS)
      TDST  = Z00*TR(IZ0,IR0,IGPOS) + Z01*TR(IZ0,IR1,IGPOS)
     &      + Z10*TR(IZ1,IR0,IGPOS) + Z11*TR(IZ1,IR1,IGPOS)
      XACT = 0.

C   Source function for dust (Kirchhoff)
C
      CALL PLANKF(TDST,FQ,BF,NF)

C   Correct emissivity for background contribution..
C
      DO L=1,NF
         EXT1(L)=gg0
         SNU (L)=gg0
         BF  (L)=BF(L)-BBOUT(L)
      ENDDO
      
C  Dust sublimation:
C
      DO ID=1,NDST
         ZV(ID)=MAX(gg0,(gg1-TDST/TSUB(ID))*20.)+1.E-10
         ZV(ID)=MIN(gg1,ZV(ID))
      ENDDO
      
      DO ID=1,NDST
         DEN=(GM(ID)+GM(1)*(gg1-ZV(1))*ERH(ID))*AB(ID)*RHOACT
         DO L=1,NF
C  Contributions of the individual dust types to the extinction:
            DEXT   =ZV(ID)*FEE(L,ID)*DEN
            EXT1(L)=EXT1(L)+DEXT
C  Contributions of the individual dust types to the emission:
            SNU (L)=SNU(L)+BF(L)*DEXT*(gg1-FEA(L,ID))
         ENDDO
         if(NFz.gt.0) then
            DO K=1,NFz
            SMzz = Z00*SMz(IZ0,IR0,IGPOS,K) + Z01*SMz(IZ0,IR1,IGPOS,K)
     &           + Z10*SMz(IZ1,IR0,IGPOS,K) + Z11*SMz(IZ1,IR1,IGPOS,K)
               L=LzMAP(K)
               DEXT   =ZV(ID)*FEE(L,ID)*DEN
C  Contributions of the individual dust types to the scattering:
               SNU(L)=SNU(L) + SMzz*DEXT*FEA(L,ID)
            ENDDO
         endif
      ENDDO
C
      DO L=1,NF
         SNU(L)=SNU(L)/(EXT1(L)+ggK)
         EXT(L)=EXT1(L)
      ENDDO

C  If non-negligible ionization then figure out its contribution

      if(xact.gt..001) then

      TACT15=TACT**1.5

C   obtain source funct.(Kirchhoff)
C
      CALL PLANKF(TACT,FQ,BF,NF)
C
C   Absorptionskoeffizient f"ur Frei-frei-Strahlung nach Spitzer:
C
      DO 30 L=1,NF
         ENP=((RHOACT/HMAS)*XACT)**2
         GFF=max(gg1,.5513*(LOG(TACT15/FQ(L))-6.4))
         EXT(L)=1.9798E-23*GFF*ENP/(TACT15*FQ(L)**2)
 30   CONTINUE
C
C   Combination of both emission from both dust and ionized gas:
C
      DO 70 L=1,NF
         SNU(L)=EXT(L)*BF(L)+EXT1(L)*SNU(L)
         EXT(L)=EXT(L)+EXT1(L)
         SNU(L)=SNU(L)/EXT(L)
 70   CONTINUE
         
      endif
    
 98   CONTINUE
C
C   Determine the next point for evaluation:
C
      CALL POSBEST(U2,IENDE)

      RETURN
      END

      SUBROUTINE POSBEST(U2,IENDE)
C This routine gives the position of the next evaluation point in
C the 2-D nested grids. The necessary variables for the subroutine
C SOURCE are in the common block POSITION.
C IGPOS: present grid
C ZPOS:  z-coordinate of evaluation site in "smallest grid units"
C RPOS:  r-coordinate "       "      "    "    "        "    "
C DUPOS: distance between old and new evaluation site
C RRPOS: current value for Rpos**2 + Zpos**2
C RRzone: central zone of special consideration when r**2+z**2 < RRzone
C DUzone: value for DUPOS within central zone
C IZPOS  Z-index for present grid
C IRPOS: R-index for present grid 
C U2:    new position on LOS
C IENDE: flag for end of LOS-integration (boundary of grid reached)

      INCLUDE 'tracedefs.h'

      data gg0,gg1/0.,1./
      data icount/-1/
      data ifirst/0/

      save

      if(ifirst.eq.0) then
        ifirst=1
        SINt =A(8)
        COSt =A(9)
        SINIt=A(10)
        COSIt=A(11)
        RRz=(RRzone/SCALZ(MGG))**2
        DUz=DUzone/SCALZ(MGG)
      endif

      if(icount.ge.0) then
        icount=icount+1
        if(icount.ge.500) stop 'icount'
        write(6,*) 'SIN,COS,SINI,COSI:',SINt,COSt,SINIt,COSIt
        write(6,666) UPOS,RPOS,ZPOS,IRPOS,IZPOS,IGPOS
 666    format('U,R,Z,IR,IZ,IG:',1p3e10.2,3i5)
      endif

      if(IENDE.ne.-111) THEN
C  Get next point, if UPOS has to be updated
        DU=.49*2**(MGG-IGPOS)
C  Special treatment (small steps) necessary when close to origin.
        if(RPOS**2 + ZPOS**2.lt.RRz) DU=DUz
        UNEXT=UPOS+DU
        DUPOS=(UNEXT-UPOS)*SCALZ(MGG)
        UPOS=UNEXT
      endif

      ZPOS=Z0+UPOS*COSt
      RPOS=sqrt(X0**2+(Y0-SINt*UPOS)**2)
      U2=UPOS*SCALZ(MGG)

      ZC=ABS( ZPOS )
      RAND=float(MZ-4)
      DO IG=MGG,1,-1
        IF(RPOS.LT.RAND .AND. ZC.LT.RAND) GOTO 6
        RPOS=RPOS*.5
        ZC=ZC*.5
      ENDDO
      RPOS=2.*RPOS
      ZC  =2.*ZC
      IG=1
C  Check if outside maximum radius
      IENDE=1
      IF(RPOS.gt.RAND+1. .or. ZC.gt.RAND+1.) GOTO 988
  6   IGPOS=IG 
      IENDE=0
      IRPOS=RPOS+.5
      IZPOS=ZC+.5
      IR0 = max(2,IRPOS+1)
      IR1 = IRPOS+2
      IZ0 = max(2,IZPOS+1)
      IZ1 = IZPOS+2
      WS1 = RPOS - float(IRPOS) +.5
      WS0 = 1. - WS1
      WT1 = ZC - float(IZPOS) +.5
      WT0 = 1. - WT1
      Z00 = WT0*WS0
      Z01 = WT0*WS1
      Z10 = WT1*WS0
      Z11 = WT1*WS1

 988  CONTINUE
      RETURN
      END

C ************************************************************************


      SUBROUTINE INITLIN(XP0,YP0,ZP0,U2,IENDE)

      INCLUDE 'tracedefs.h'

      DUPOS=0.

C Convert the LOS parameters X0,Y0 and Z0 to "smallest grid units"
      X0=XP0/SCALZ(MGG)
      Y0=YP0/SCALZ(MGG)
      Z0=ZP0/SCALZ(MGG)
        
C  Maximum values for abs(Upos), R and Z:
      RAND=(MR-3.5)*2**(MGG-1)
      if(abs(X0).gt.RAND .or.
     *   abs(Y0).gt.RAND .or.
     *   abs(Z0).gt.RAND) goto 9999

      RRAND=RAND**2
      SINt =A(8)
      COSt =A(9)
      SINIt=A(10)
      COSIt=A(11)

C  Determine whether to start from the outer boundary or from
C  the plane  Zpos=0.

      if(IENDE.eq.0) then
C  Here we start from Zpos=0 
        UPOS=-Z0*COSIt
        YTEST=Y0-UPOS*SINt
c       write(6,*) 'Ytest=',YTEST
        if(X0**2+YTEST**2.gt.RRAND) then
          YTEST=SIGN(sqrt(RRAND-X0**2),YTEST)
c         write(6,*) 'Ytest=',YTEST
          UPOS=-abs((Y0-YTEST)*SINIt)
        endif
      else
C  Here we start from the outer boundary Upos = -Rand
C  First, estimate the values...
        B=Z0*COSt-Y0*SINt
        C=RRAND-X0**2-Y0**2-Z0**2
        if(C.lt.0.) goto 9999
        UPOS=-B-sqrt(B**2+C)
      endif

      IENDE=-111
      CALL POSBEST(U2,IENDE)
      RETURN

C Set flag IENDE, if LOS is outside the outer grid.
 9999 IENDE=999
      RETURN
      END
C
      SUBROUTINE STATEM(E,RHO,X,PP,YY,TT,EE)
C inputs:
C E : specific energy
C RHO : Density
C X : ionization degree of hydrogen
C
C outputs:
C PP : gas pressure
C TT : gas temperature
C workspace
C YY : mass fraction of atomic hydrogen
C EE : specific energy as returned from ESV

      INCLUDE 'tracedefs.h'
      PARAMETER ( FF=BOLTZ/HMAS , FF32=1.5*FF , FF54=1.25*FF ,
     &    FF12=0.5*FF , FF14=0.25*FF , EBND12=0.5*EBNDH2/HMAS )
      PARAMETER (MSE=300,MSR=80)
      DIMENSION E(MZZ,MRR),RHO(MZZ,MRR),PP(MZZ,MRR),
     $     YY(MZZ,MRR),TT(MZZ,MRR),EE(MZZ,MRR),X(MZZ,MRR)

      COMMON /ST/TL(MSE,MSR)
      SAVE
      DATA RLMIN,RLMAX,ELMIN,ELMAX/-20.,-5.,10.9,12.8/
      DATA IFIRST/0/
C
C  INITIALIZE STATE IF FIRST CALL
C
      IF(IFIRST.NE.1) THEN
        IFIRST=1
        UM=10.
        UM=LOG(UM)
        DRL=(RLMAX-RLMIN)/FLOAT(MSR-1)
        DEL=(ELMAX-ELMIN)/FLOAT(MSE-1)
        EMIN=EXP(UM*ELMIN)
        EMAX=EXP(UM*ELMAX)
        TS=EMIN/FF54
        DO 4 IR=1,MSR
          RHL=RLMIN+DRL*FLOAT(IR-1)
          RH=EXP(UM*RHL)

          T2=TS
          DO 3 IE=1,MSE
            EL=ELMIN+DEL*FLOAT(IE-1)
            EN=EXP(UM*EL)

            DO 1 IT=1,50

              CALL EST(RH,T2,EPS,YYA,EPT)

C              WRITE(6,766)  IT,RH,T2,EPS,YYA,EPT
C 766          FORMAT(' it=',i3,' rho,t,eps,yya,ept=',1p,5e10.2)
              DT=(EN-EPS)/EPT
              T2=T2+DT
              IF(ABS(DT/T2).LT.1.E-6) GOTO 2
 1          CONTINUE
            WRITE(6,435) IR,IE,RHL,RH,EL,EN
 435        FORMAT(' NO CONVERGENCE IN STATE SETUP: IR,IE=',2I4,
     &      '  RHL,RH=',0PF8.3,1PE12.4,' EL,E=',0PF8.3,1PE12.4)
            STOP
 2          TL(IE,IR)=LOG10(T2)
 3        CONTINUE
 4      CONTINUE
      ENDIF
C  END OF STATE INITIALIZATION
C
      DO 40 IR=1,MR
      DO 40 IZ=1,MZ
        IF(E(IZ,IR).GT.EMAX) THEN
          TT(IZ,IR)=(E(IZ,IR)-EBND12)/FF32
        ELSE
          IF(E(IZ,IR).LT.EMIN) THEN
            TT(IZ,IR)=E(IZ,IR)/FF54
          ELSE
            RL=(LOG10(RHO(IZ,IR))-RLMIN)/DRL
            J1=RL+1.
            J1=MAX(J1+1,1)
            J1=MIN(J1,MSR-1)
            J2=J1+1
            RL=RL-FLOAT(J1-1)
            EL=(LOG10(E  (IZ,IR))-ELMIN)/DEL
            K1=EL+1.
            K1=MAX(K1+1,1)
            K1=MIN(K1,MSE-1)
            K2=K1+1
            EL=EL-FLOAT(K1-1)
            TT(IZ,IR)=EXP(UM*((1.-RL)*((1.-EL)*TL(K1,J1)+EL*TL(K2,J1))
     &                     + RL   *((1.-EL)*TL(K1,J2)+EL*TL(K2,J2)) ))
          ENDIF
        ENDIF
 40   CONTINUE
C
* ORIGINAL VERSION FOR ONLY NEUTRAL OR MOLECULAR GAS
*     CALL ESV(RHO,TT,EE,YY,PP)
C      dtmax=0.
*     DO 50 IR=2,MR-1
*     DO 50 IZ=2,MZ-1
C add a final Newton-Raphson step for temperature
*       DT=(E(IZ,IR)-EE(IZ,IR))/PP(IZ,IR)
*       TT(IZ,IR)=TT(IZ,IR)+DT
C calculate gas pressure
*       PP(IZ,IR)=FF12*RHO(IZ,IR)*TT(IZ,IR)*(1.+YY(IZ,IR))
C      if(abs(dt/tt(iz,ir)) .lt. abs(dtmax)) goto 50
C        irs=ir
C        izs=iz
C        dtmax=dt/tt(iz,ir)
*50   CONTINUE
C        write(6,333) izs,irs,dtmax,dtmax*tt(izs,irs),tt(izs,irs)
C 333    format(' max. korr. at (',i3,',',i3,') dtmax,dt,tt=',1p,3e12.4)
*
* NEW VERSION INCLUDING IONIZATION (A.W. JUNE 1992)
C
      Z0=0.0
      TMIN=3.0
      CALL ESV(RHO,TT,EE,YY,PP)
      DO 59 IR=1,MR
      DO 59 IZ=1,MZ
         DT=(E(IZ,IR)-EE(IZ,IR))/PP(IZ,IR)
         TT(IZ,IR)=TT(IZ,IR)+DT
   59 CONTINUE
      CALL ESV(RHO,TT,EE,YY,PP)
      DO 60 IR=1,MR
      DO 60 IZ=1,MZ
         FAC=MIN( X(IZ,IR) , XL )
         FAC = ( XLA*FAC + XLB )*FAC**2 + 1.0
         YY(IZ,IR)= 1.0 - (1.0-YY(IZ,IR))*FAC
         TT(IZ,IR)= (E(IZ,IR)-YY(IZ,IR)*EBND12) /
     &     ( FF54 + YY(IZ,IR)*FF14 + X(IZ,IR)*FF32 )
         TT(IZ,IR)= MAX( TT(IZ,IR) , TMIN )
         PP(IZ,IR)= FF12*RHO(IZ,IR)*TT(IZ,IR)*
     &             ( 1. + YY(IZ,IR) + 2.*X(IZ,IR) )
   60 CONTINUE
C
      RETURN
      END
C
C ****************************************************************
C
      SUBROUTINE EST(RH,TM,EPS,YYA,EPT)
      INCLUDE 'tracedefs.h'

      PARAMETER ( FF=BOLTZ/HMAS , FF32=1.5*FF , FF54=1.25*FF ,
     &    FF12=0.5*FF , FF14=0.25*FF , EBND12=0.5*EBNDH2/HMAS )
      SAVE
C INPUT:
C  RH= DENSITY
C  TM= GAS TEMPERATURE
C OUTPUT
C   YYA=MASS FRACTION OF ATOMIC HYDROGEN
C  EPS=specific energy
C  EPT=partial derivative dEPS/dTM
C
        IF(TM.LT.700.) THEN
          YYA=0.
          EPS=FF54*TM
          EPT=FF54
        ELSE
          AH=RH/2.11*EXP(52490./TM-20.*RH)
          B=SQRT(1.+4.*AH)
C YY=mass fraction of atomic hydrogen
          YYA=2./(1.+B)
C YYT= partial derivative dYY/dTM
          YYT=52490. * AH/B * (YYA/TM)**2
C EPS=specific energy
          EPS=FF54*TM + YYA * (FF14*TM+EBND12)
C EPT= partial derivative dEPS/dTM
          EPT=FF54+FF14*YYA + YYT * (FF14*TM+EBND12)
        ENDIF
      RETURN
      END
C
C **********************************************
C
      SUBROUTINE ESV(RH,TM,EPS,YY,EPT)
C     vector version of EST
      INCLUDE 'tracedefs.h'

      PARAMETER ( FF=BOLTZ/HMAS , FF32=1.5*FF , FF54=1.25*FF ,
     &    FF12=0.5*FF , FF14=0.25*FF , EBND12=0.5*EBNDH2/HMAS )
      DIMENSION RH(MZZ,MRR),TM(MZZ,MRR),EPS(MZZ,MRR),
     $     YY(MZZ,MRR),EPT(MZZ,MRR)
      SAVE
C
      DO 20 IR=1,MRR
      DO 20 IZ=1,MZZ
        IF(TM(IZ,IR).LT.700.) THEN
          YY(IZ,IR) = 0.
          EPS(IZ,IR)= FF54*TM(IZ,IR)
          EPT(IZ,IR)= FF54
        ELSE
          AH = RH(IZ,IR)/2.11*EXP(52490./TM(IZ,IR)-20.*RH(IZ,IR))
          B = SQRT(1.+4.*AH)
          YY(IZ,IR) = 2./(1.+B)
          YYT       = 52490. * AH/B * (YY(IZ,IR)/TM(IZ,IR))**2
          B         = FF14*TM(IZ,IR) + EBND12
          EPS(IZ,IR)= FF54*TM(IZ,IR) + YY(IZ,IR)*B
          EPT(IZ,IR)= FF54+FF14*YY(IZ,IR) + YYT*B
        ENDIF


  20  CONTINUE
      RETURN
      END
