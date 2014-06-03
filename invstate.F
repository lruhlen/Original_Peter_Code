      subroutine invstate(im,Pr,temp,X,Y,Crad,rho,cP,alpha,Rbeta1,delta)
C+
C  The subroutine invstate determines the density as a function of
C  pressure and temperature.  Other thermodynamic quantities are
C  also calculated.  This current version includes the effects of
C  degeneracy, including partial relativistic degeneracy; it does
C  not include partial ionization of the degenerate gas.
C
C  Author:  Harold W. Yorke      10-SEP-02               (JPL / CALTECH)
C-
      include 'parm.h'
      real*4 Plog,Tlog,XX,YY,Ymax,EOS(5)
      real*8 XCA,H1,A
      COMMON/ABUND/XCA(14),H1(14),A(14)
      data iorder/5/
      DATA Arad3/2.52197145D-15/
      DATA negative/0/
      DATA ifirst/0/
      save
C
      if(ifirst.eq.0) then
        ifirst=1
        Ymax=H1(1)+H1(2)
      endif
      Prad=Crad*Arad3*Temp**4
      Pgas=Pr-Prad
      if(Pgas.lt.0.d0) then
        write(6,'(a,i3,a)') 'INVSTATE',im,': Gas pressure is negative'
        write(6,'(i4,1x,a,1p,4e12.4)')                                  &
     &     negative,'Pr,temp,Prad,Pgas:',Pr,temp,Prad,Pgas
        negative=negative+1
        if(negative.gt.10) stop 'invstate'
        Pgas = 0.01*Pr
        Prad = 0.99*Pr
      endif
      Rbeta=Pgas/Pr
      Rbeta1=Prad/Pr
      i1=0
c
      Plog=log10(Pgas)
      Tlog=log10(Temp)
      XX=X
      YY=Y
      call eostab(Plog,Tlog,XX,YY,Ymax,EOS,iorder,ierror)
      if(ierror.ne.0) then
c       write(6,201) Plog,Tlog,EOS,iorder,ierror
 201    format(7f10.5,2i2)
c       stop 'invstate'
      endif
      rho     = 10.**eos(2)
      alphgas = eos(3)
      deltgas = eos(4)
      cVgas   = 10.**(eos(5)+6.)
      cPgas   = cVgas+deltgas**2*Pgas/(alphgas*temp*rho)
c     atg     = deltgas*Pgas/(rho*temp*cPgas)
c
c CHANGE FINAL RESULTS TO ACCOUNT FOR RADIATION PRESSURE
c
      alpha = alphgas/Rbeta
      delta = deltgas + 4.d0*alpha*Rbeta1
      cV    = cVgas + 12.d0*Prad/(temp*rho)
      cP    = cV + delta**2*Pr/(alpha*temp*rho)
c
      return
      end
      subroutine eostab(Plog,Tlog,X,Y,Ymax,EOS,iorder,ierror)
      implicit real*4(a-h,o-z) , integer (i-n)
      parameter (NZ=5,NX=5,MT=139,NTX=(NX+NZ)*MT,MPT=187165)
      dimension IPTR(NTX+1),EOSA(MPT,5)
      dimension EOS(5)
      dimension Pmin(4),DP1(4),DP2(4),XX(4),IP(4),IP1(4),IP2(4),NP(4)
      data ifirst/0/
      save

      ierror=0
      if(ifirst.eq.0) then
        ifirst=1
        open(7,file='eospointer',form='formatted',status='old')
        read(7,*) IPTR,kx,kt,tmin,tmax,dP
        close(7)
        open(7,file='eostable',form='formatted',status='old')
c
        if(kx.ne.NX .or. kt.ne.MT) then
          write(6,*) 'Mismatch of dimensions in eostab and eostable',
     &     ' NX=',NX,kx,' MT=',MT,kt
          stop 'EOSTAB'
        endif

        dx  = 1./float(NX)
        dz  = 1./float(NZ)
        Tlogmin=log10(Tmin)
        Tlogmax=log10(Tmax)
        dt  = (Tlogmax-Tlogmin)/float(MT-1)

        ITX = 0
        do ix=1,NX+NZ
        do it=1,MT
          ITX=ITX+1
          read (7,201) (EOSA(IPTR(ITX),i),i=1,5)
 201      format(4x,1p,E12.4,7x,f9.6,5x,f4.1,3x,f9.6,3x,f9.6)
          IPmin=IPTR(ITX)+1
          IPmax=IPTR(ITX+1)-1
          do j=IPmin,IPmax
            read(7,203) (EOSA(j,i),i=1,5)
 203        format(f10.6,f11.6,3f9.6)
          enddo
        enddo
        enddo
      close(7)
      endif

      if(X.lt.1.E-37 .and. Ymax-Y.gt.1.e-4) then
        ZZ=1.-Y
        DX2=ZZ/dz
        IX1=max(DX2,0.)+4.
        IX2=min(IX1+2,NX+NZ)
        IX1=IX2-1
        if(IX1.gt.NX) then
          DX2=DX2-float(IX1-5)
        else
          IX1=1
          DX2=max(0.,(ZZ-1.+Ymax)/(dz-1.+Ymax))
        endif
        DX1=1.-DX2
      else
        if(X.lt.0.0 .or. X.gt.0.8) ierror=1
        DX2=X/dx
        IX1=max(DX2,0.)
        IX2=min(IX1+2,NX)
        IX1=IX2-1
        DX2=DX2-float(IX1-1)
        DX1=1.-DX2
      endif
c     write(6,202) 'X:',ix1,ix2,dx1,dx2
 202  format(a,2i10,1p,2E12.4)

      if(Tlog.lt.Tlogmin .or. Tlog.gt.Tlogmax) ierror=ierror+2
      DT2=(Tlog-Tlogmin)/dt
      IT1=max(DT2,0.)
      IT2=min(IT1+2,MT)
      IT1=IT2-1
      DT2=DT2-float(IT1-1)
      DT1=1.-DT2
c     write(6,202) 'T:',it1,it2,dt1,dt2

      IP(1)=(IX1-1)*MT+IT1
      IP(2)=(IX2-1)*MT+IT1
      IP(3)=(IX1-1)*MT+IT2
      IP(4)=(IX2-1)*MT+IT2
c     write(6,*) (IPTR(IP(j)),j=1,4)
      irror=0
      do j=1,4
        Pmin(j) = EOSA(IPTR(IP(j))+1,1)
        NP(j)   = IPTR(IP(j)+1)-IPTR(IP(j))-1
        DP2(j)  = (Plog-Pmin(j))/DP
        IP1(j)  = max(DP2(j),0.)
        IP2(j)  = min(IP1(j)+2,NP(j))
        IP1(j)  = IP2(j)-1
        DP2(j)  = DP2(j)-float(IP1(j)-1)
        DP1(j)  = 1.-DP2(j)
        if(DP1(j).lt.0.0 .or. DP2(j).lt.0.0) irror=4
      enddo
      ierror=ierror+irror
c     do j=1,4
c       write(6,202) 'P:',iP1(j),iP2(j),dP1(j),dP2(j)
c     enddo

      do i=1,iorder
        do j=1,4
        XX(j)=DP1(j)*EOSA(IPTR(IP(j))+IP1(j),i) + 
     &        DP2(j)*EOSA(IPTR(IP(j))+IP2(j),i)
        enddo
        if(i.le.2) then
          DS1=DT1
          DS2=DT2
        else
          DS1=min(1.,DT1)
          DS2=1.-DS1
        endif
        EOS(i)=DX1*(DS1*XX(1)+DS2*XX(3)) + DX2*(DS1*XX(2)+DS2*XX(4))

      enddo

      return
      end
