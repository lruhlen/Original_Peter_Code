      SUBROUTINE HENYEY(MODE)
      include 'parm.h'
      include 'var.h'
      logical LCONV
      DIMENSION D(MH,MH+1),CA(MH),DX(MH)
      EQUIVALENCE (HA(1,MH+1),D(1,1)) , (D(1,MH+1),CA(1))
      save
C
      if(abs(ITMIN).gt.10 .or. mode.eq.999) then
        cf=1./max(10,ITMIN)
        cfmax=.6
        cfexp=1.2
        if(ITMIN.gt.10) cfexp=1.1
      else
c       cf=.1
        cf=1.
        cfmax=1.0
        cfexp=1.25
      endif
      ID=N
      MODE=0
      IRIT=0
      JTMAX=abs(ITMAX)
      IF(ITMIN.GT.0 .and. ITMAX.gt.0) IRIT=JTMAX/2
      ITER=0
  300 CONTINUE
      if(iter.eq.ITMAX/2) cf=.1
      cf=min(cfmax,cf*cfexp)
      ITER=ITER+1
      IF(ITER.GT.JTMAX) GOTO 998
        J2=1
        J3=2
        CALL GI(J2)
        IF(ITMIN.EQ.-ITER) then
          call gid(j2)
          WRITE(6,101) 'J2=',J2,' Henyey matrices'
          WRITE(6,888) 'G:',(G(I),I=1,NG)
          do j=1,ng
            WRITE(6,888) 'C:',(HC(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'D:',(HD(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'E:',(HE(I,J),I=1,NG)
          enddo
        endif
        do I=1,MH
          do J=1,MH
            HA(I,J)=HD(I,J)
            D(I,J)=-HE(I,J)
          enddo
          D(I,MH+1)=-G(I)
          DG(I)=0.
          JG(I)=0
        enddo
        CALL GIRL (HA,MH,MH+1)
        IW=1
        do I=1,MH
          do J=1,MH+1
            HX(I,J)=D(I,J)
            HW(I,J,IW)=D(I,J)
          enddo
        enddo
  302 CONTINUE
        J1=J2
        J2=J3
        J3=J3+1
        CALL GI(J2)
        IF(ITMIN.EQ.-ITER) then
          call gid(j2)
          WRITE(6,101) 'J2=',J2,' Henyey matrices'
          WRITE(6,888) 'G:',(G(I),I=1,NG)
          do j=1,ng
            WRITE(6,888) 'C:',(HC(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'D:',(HD(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'E:',(HE(I,J),I=1,NG)
          enddo
        endif
        do I=1,MH
          do J=1,MH
            FAC=HC(I,1)*HX(1,J)
            do K=2,MH
              FAC=FAC+HC(I,K)*HX(K,J)
            enddo
            HA(I,J)=HD(I,J)+FAC
            D(I,J)=-HE(I,J)
          enddo
          FAC=HC(I,1)*HX(1,MH+1)
          do K=2,MH
            FAC=FAC+HC(I,K)*HX(K,MH+1)
          enddo
          D(I,MH+1)=-G(I)-FAC
          if(ABS(DG(I)).lt.ABS(G(I))) then
            DG(I)=G(I)
            JG(I)=J2
          endif
        enddo
        CALL GIRL (HA,MH,MH+1)
        IW=IW+1
        do J=1,MH+1
          do I=1,MH
            HX(I,J)=D(I,J)
            HW(I,J,IW)=D(I,J)
          enddo
        enddo
      IF(J2.GE.N-1) GOTO 303
      IF(IW.LT.ID) GOTO 302
        IW=0
      GOTO 302
  303 IF(ITER.GT.IRIT .and. mod(MODEL,JRIT).eq.0)
     *  WRITE(6,101) 'ITER:',ITER,' LARGEST GI:',(JG(I),DG(I),I=1,MH)
        J1=J2
        J2=J3
        J3=J2+1
        CALL GI(J2)
        IF(ITMIN.EQ.-ITER) then
          call gid(j2)
          WRITE(6,101) 'J2=',J2,' Henyey matrices'
          WRITE(6,888) 'G:',(G(I),I=1,NG)
          do j=1,ng
            WRITE(6,888) 'C:',(HC(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'D:',(HD(I,J),I=1,NG)
          enddo
          do j=1,ng
            WRITE(6,888) 'E:',(HE(I,J),I=1,NG)
          enddo
        endif
        do I=1,MH
          do J=1,MH
            FAC=HC(I,1)*HX(1,J)
            do K=2,MH
              FAC=FAC+HC(I,K)*HX(K,J)
            enddo
            D(I,J)=HD(I,J)+FAC
          enddo
          FAC=HC(I,1)*HX(1,MH+1)
          do K=2,MH
            FAC=FAC+HC(I,K)*HX(K,MH+1)
          enddo
          CA(I)=-G(I)-FAC
          G(I)=DG(I)
          DG(I)=0.
          JG(I)=0
        enddo
        CALL GIRL (D,MH,1)
  304 CONTINUE
        do I=1,MH
          E(I)=CA(I)
          AFAC=E(I)*cf
          IF(SMAX(I).gt.0.) then
            cm39 = ABS(AFAC/X(j2,I))
            if (cm39 .gt. SMAX(I)) then
              AFAC = AFAC/ABS(AFAC) * SMAX(I)*X(J2,I)
            end if 
            X(J2,I)=X(J2,I)+AFAC   
            DX(I)=2.*AFAC/(X(J2,I)+VX(J2,I)+SMIN(I))
            IF(X(J2,I).lt.SMIN(I)) goto 997
          else
            AFAC=MIN(AFAC,-SMAX(I))
            AFAC=MAX(AFAC, SMAX(I))
            X(J2,I)=X(J2,I)+AFAC
            DX(I)=2.*AFAC/(ABS(X(J2,I))+ABS(VX(J2,I))+SMIN(I))
          endif
          IF(ABS(DG(I)).lt.ABS(DX(I))) then
            DG(I)=DX(I)
            JG(I)=J2
          endif
        enddo
        IF(ITMIN.LT.0) then
          WRITE(6,101) 'J2=',J2,' CORRECTIONS:',(I,E(I),I=1,MH)
        ENDIF
      IF(J2.EQ.1) GOTO 305
        do I=1,MH
          CA(I)=HW(I,MH+1,IW)
          do J=1,MH
            CA(I)=CA(I)+HW(I,J,IW)*E(J)
          enddo
        enddo
        J2=J1
        J1=J1-1
        IW=IW-1
      IF(IW.GT.0 .OR. J2.EQ.1) GOTO 304
        IW=ID
      GOTO 304
  305 CONTINUE
        IF(ITER.GT.IRIT .and. mod(MODEL,JRIT).eq.0) then
          WRITE(6,101) 'ITERATION',ITER,' CORRECTIONS:',
     *    (JG(I),DG(I),I=1,MH)
          WRITE(6,101) 'ITERATION',ITER,'  NEW VALUES:',
     *    (JG(I),X(J2,i),i=1,MH)
        ENDIF
        LCONV=.FALSE.
        do I=1,MH
          LCONV=LCONV .or. ABS(DG(I)) .gt. EPS(I)
        enddo
      IF(ITER.LT.ITMIN .OR. LCONV .OR. abs(cf-cfmax).gt.0.01d0) GOTO 300
        do I=1,MH
          DG(I)=0.
          JG(I)=0
        enddo
        do J=2,N
          do I=1,MH
            E(I)=2.*(X(J,I)-VX(J,I))/(ABS(X(J,I))+ABS(VX(J,I))+SMIN(I))
            if(ABS(DG(I)).lt.ABS(E(I))) then
              DG(I)=E(I)
              JG(I)=J
            endif
          enddo
        enddo
        if(mod(MODEL,JRIT).eq.0) then
          WRITE(6,101) 'MODEL ',MODEL,' CHANGES:',(JG(I),DG(I),I=1,MH)
          WRITE(6,101) 'ITER=',ITER,' GI:',(I,G(I),I=1,MH)
        endif
      GOTO 999
 997    MODE=I
        WRITE(6,103) J2,I,X(J2,I)
      GOTO 999
 998    MODE=999
        WRITE(6,101) '--- CONVERGENCE NOT POSSIBLE: ',MODEL,' ',
     *               (JG(I),DG(I),I=1,MH),ITER
 999  CONTINUE
      CHANGE=0.
      DO I=1,MH
        CHANGE=CHANGE+ABS(DG(I))
      ENDDO
      RETURN
  101 FORMAT(1X,A,i5,A,1P,7(I5,E10.2),/,7x,8(I5,E10.2))
  103 FORMAT(/6X,'VARIABLE BELOW LIMIT:      X(',I4,',',I1,') = ',
     1       1P,E12.5)
  888 FORMAT(8x,a,1P,1x,10E10.2,10(/1X,10E10.2))
      END
      SUBROUTINE GIRL(A,N,M)
      include 'parm.h'
      DIMENSION A(1)
      SAVE
      NPM = N+M
      DO J = 1,N
        NJ = (J-1)*N
        JJ = NJ + J
        J1 = J + 1
        AMAX =ABS (A(JJ))
        JM = J
        IF( J1 .le. N ) THEN 
          DO I = J1,N
            IJ = NJ + I
            IF(ABS (A(IJ)) .gt. AMAX ) THEN
              AMAX =ABS (A(IJ))
              JM = I
            ENDIF
          ENDDO 
          IF( JM .ne. J ) THEN
            I1 = JM + NJ
            I2 = JJ
            DO I = J,NPM
              ZWI = A(I1)
              A(I1) = A(I2)
              A(I2) = ZWI
              I1 = I1 + N
              I2 = I2 + N
            ENDDO
          ENDIF
        ENDIF
        IF( A(JJ) .eq. 0.) GOTO 40
        DO I = 1,N
          IF( I .ne. J ) THEN
            IJ = NJ + I
            IK = NJ + I
            JK = JJ
            FAKTOR = - A(IJ)/A(JJ)
            DO K = J1,NPM
              JK = JK + N
              IK = IK + N
              A(IK) = A(IK) + FAKTOR * A(JK)
            ENDDO
          ENDIF
        ENDDO
        JK = JJ
        FAKTOR = 1./A(JJ)
        DO K = J1,NPM
          JK = JK + N
          A(JK) = A(JK) * FAKTOR
        ENDDO
      ENDDO
      RETURN
   40 CONTINUE
      WRITE(6,100)
  100 FORMAT(' ERROR EXIT IN SUBROUTINE GIRL')
      STOP 'girl'
       END
