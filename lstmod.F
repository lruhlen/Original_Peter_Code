      PROGRAM LSTMOD
C+
C  THE PROGRAM LSTMOD PRINTS OUT MODEL INFORMATION
C
C  AUTHOR: H.W. YORKE   27-OCT-02               (JPL / CALTECH)
C-
      include 'parm.h'
      include 'var.h'
      character*10 SFILE
      character*1 OFILE(10)
      equivalence (SFILE,OFILE)
      DATA MODWRT/-1/
      save
C
      Crad=1.d0
      do i=1,999999
        READ (9,END=999) MODEL,N,TIME,dTIME,TAUJK,                      &
     &     ((x(j,k),k=1,4),dM(j),                                       &
     &     hydrogen(j),helium3(j),helium4(j),deuterium(j),              &
     &     carbon(j),trogen(j),oxygen(j),                               &
     &     j=1,N)
        zM(1)=dM(1)
        do j=2,N
          zM(j)=zM(j-1)+dM(j)
        enddo
        write(SFILE,200) i
 200    format('star',i6)
        do j=5,10
          if(OFILE(j).eq.' ') OFILE(j)='0'
        enddo
        open(unit=8,file=SFILE,form='formatted',status='new') 
        do j=1,N
          call invstate(i1,x(j,1),x(j,4),hydrogen(j),helium4(j),Crad,
     &         RHO,cP,alpha,delta)
          if(j.ne.N) then
            energ=(X(j+1,3)-X(j,3))/dM(j+1)
          else
            energ=0.
          endif
          write (8,202) J,ICONVECT(j),dM(J),zM(J),
     *                  (X(j,k),k=1,MH),rho,energ,
     x                  hydrogen(j),helium3(j),helium4(j),deuterium(j),
     *                  carbon(j),oxygen(j)
202       format(I5,I2,1p,14E12.4)
        enddo
      enddo
 999  stop
      end
