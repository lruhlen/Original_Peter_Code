C=======================================================================
C
C
C                        INCLUDE FILE treedefs.h                       
C
C
C=======================================================================
C
C
C     Parameter declarations, allocation of array storage, common
C     block definitions.
C
C
C=======================================================================

        CHARACTER*9 boundary
        CHARACTER*50 headline
        INTEGER root,subp,nbodies,incells,nttot,ntmin,ntmax,ntavg,
     &          nsteps,noutbod,noutlog,ndim,nsubcell,nbodsmax,ncells,
     &          nbodcell,nbods1,incellsg,nsat,intzero,intone,inttwo,
     &          intthree,intfour,nsvolume,nstot,nsmin,nsmax,nsavg
        LOGICAL usequad,usesph,fixthalo,selfgrav,fixtriax,variabls,
     &          dgravsum,usebh,inpteps,restart,rigidbar
        REAL mass,tol,tol2,eps,rsize,rmin,phi,pos,vel,acc,quad,tnow,
     &       tpos,dtime,dtime2,tiny,one,two,four,cputime0,cputime1,etot,
     &       cputime,zero,zero02,cellsize,minusone,ttree,mtot,ektot,
     &       eptot,three,log2,mstar,mgas,snheat,minustwo,vhalo,gamhalo,
     &       rhalo,atriax,btriax,ctriax,mastriax,epsvect,onehalf,
     &       massbh,phiext,tstartbh,etabh,nsvtol,posbh,velbh,amvecbh,
     &       cmvelbh,cmposbh,esofttot,patspeed,barmass,barsize,
     &       tgrowbar,tlivebar,tdiebar,enerror

        PARAMETER(ndim=3,nsubcell=2**ndim)
        PARAMETER(nbodsmax=65536,ncells=40000)
        PARAMETER(nbodcell=nbodsmax+ncells,nbods1=nbodsmax+1)

        COMMON/paramcom/nbodies,tol,tol2,eps,usequad,usesph,fixthalo,
     &                  selfgrav,variabls,nsvolume,nsvtol,dgravsum,
     &                  inpteps,restart
        COMMON/boundcom/boundary
        COMMON/msgcom/headline
        COMMON/cellcom/rsize,rmin(ndim),incells,incellsg
        COMMON/pointers/root,subp(nbods1:nbodcell,1:nsubcell)
        COMMON/bodycell/mass(1:nbodcell),phi(1:nbodsmax),
     &                  pos(1:nbodcell,1:ndim),cellsize(1:nbodcell),
     &                  vel(1:nbodsmax,1:ndim),acc(1:nbodsmax,1:ndim),
     &                  epsvect(1:nbodcell),phiext(1:nbodsmax)
        COMMON/quadcom/quad(nbods1:nbodcell,1:2*ndim-1)
        COMMON/forcecom/nttot,ntmin,ntmax,ntavg
        COMMON/softcom/nstot,nsmin,nsmax,nsavg
        COMMON/timecom/nsteps,noutbod,noutlog,tnow,tpos,dtime,dtime2,
     &                 ttree
        COMMON/cpucom/cputime0,cputime1,cputime
        COMMON/misccom/tiny,zero,one,two,four,minusone,zero02,three,
     &                 log2,minustwo,intone,inttwo,intthree,intfour,
     &                 intzero,onehalf
        COMMON/enrgycom/mtot,etot,ektot,eptot,mstar,mgas,snheat,
     &                  esofttot,enerror
        COMMON/halocom/vhalo,gamhalo,rhalo,nsat
        COMMON/triaxcom/fixtriax,atriax,btriax,ctriax,mastriax
        COMMON/barcom/rigidbar,patspeed,barmass,barsize,tgrowbar,
     &                tlivebar,tdiebar
        COMMON/blakhcom/usebh,massbh,tstartbh,etabh,posbh(1:ndim),
     &                  velbh(1:ndim),amvecbh(1:ndim),cmvelbh(1:ndim),
     &                  cmposbh(1:ndim)

C=======================================================================
C   Definitions specific to input/output.
C=======================================================================
        INTEGER uterm,upars,ulog,ubodsin,ubodsout,uboddump,utermfil,
     &          ireclog,ucrash,uindump
        CHARACTER*8 parsfile,logfile,ibodfile,obodfile,dumpfile,
     &              termfile,crashfil,indumpf

        PARAMETER(uterm=6,upars=10,ulog=11,ubodsin=12,ubodsout=13,
     &            uboddump=14,utermfil=15,ucrash=16,uindump=17)
        PARAMETER(parsfile='TREEPAR',logfile='TREELOG',
     &            ibodfile='TREEBI',obodfile='SNAPxxxx',
     &            dumpfile='SYSDUMP',termfile='TREEOUT',
     &            crashfil='CDUMP',indumpf='DUMPIN')

        COMMON/iocommon/ireclog

C=======================================================================
C   Definitions specific to individual particle time steps.
C=======================================================================
        INTEGER pactive,npactive,ntvector,itimestp,otimestp,inittbin,
     &          nsphact,mintstep,upbin
        LOGICAL endstep
        REAL etol,stime,tsteppos

        COMMON/actcom/npactive,nsphact,pactive(nbodsmax)
        COMMON/stepcom/itimestp(nbodsmax),otimestp(nbodsmax),inittbin,
     &                 mintstep,upbin,etol,stime,tsteppos,endstep,
     &                 ntvector
C=======================================================================
C   Definitions specific to repeated restarts:
C=======================================================================
        INTEGER nbig,deathlist
        REAL hlength
        PARAMETER (nbig=11)
        COMMON/rblock/hlength(nbig),deathlist(nbodsmax)
C=======================================================================
C   Definitions specific to SPH.
C=======================================================================
        CHARACTER*2 symmetry
        CHARACTER*4 artfvisc
        INTEGER nsphmax,ninterp,nsph,nsmooth,nstar
        LOGICAL sphinit,variablh,starform,friction,readinas,ethinit,
     &          readinet,outputet,outputas,uentropy,entinit,variablg,
     &          geogradp,sphfiltr,friclucy,inpmetal,inptform,isotherm
        REAL rho,ethermal,csound,ethold,dethdt,dethold,hsmooth,piinv,
     &       deldr2i,alpha,beta,epssph,gamma,gamma1,gamma2,courant,
     &       wsmooth,dwsmooth,phsmooth,acsmooth,veltpos,mumaxdvh,teth,
     &       ethtot,tvel,hsmdivv,epsgas,ctcutoff,hsmcurlv,consthsm,
     &       cfrict,entropy,entold,dentdt,dentold,epsfiltr,tempsphv,
     &       thirty5,nsmtol,tform,metals
     
        PARAMETER(nsphmax=32768,ninterp=30000)

        DIMENSION entropy(nsphmax),entold(nsphmax),dentdt(nsphmax),
     &            dentold(nsphmax)

        EQUIVALENCE (entropy(1),ethermal(1)),(entold(1),ethold(1)),
     &              (dentdt(1),dethdt(1)),(dentold(1),dethold(1))

        COMMON/statecom/rho(nsphmax),csound(nsphmax),ethermal(nsphmax),
     &                  tform(nbodsmax),metals(nbodsmax),inpmetal,
     &                  inptform,isotherm
        COMMON/ethcom/ethold(nsphmax),dethdt(nsphmax),dethold(nsphmax),
     &                teth,ethtot
        COMMON/tsynccom/veltpos(1:nsphmax,1:ndim),tvel(nsphmax)
        COMMON/artfvcom/artfvisc
        COMMON/symmcom/symmetry
        COMMON/visccom/alpha,beta,epssph,mumaxdvh(nsphmax),
     &                 hsmdivv(nsphmax),hsmcurlv(nsphmax)
        COMMON/sphparam/nsph,gamma,gamma1,gamma2,courant,sphinit,nstar,
     &                  variablh,epsgas,starform,friction,readinas,
     &                  ethinit,readinet,outputet,ctcutoff,outputas,
     &                  cfrict,uentropy,entinit,variablg,geogradp,
     &                  friclucy,thirty5
        COMMON/filtrcom/sphfiltr,epsfiltr,tempsphv(nsphmax,3)
        COMMON/smoocom/hsmooth(1:nsphmax),nsmooth,piinv,consthsm,
     &                 nsmtol
        COMMON/interpoc/deldr2i
        COMMON/skerncom/wsmooth(0:1+ninterp),dwsmooth(0:1+ninterp)
        COMMON/gravcom/phsmooth(0:1+ninterp),acsmooth(0:1+ninterp)

C=======================================================================
C   Definitions specific to radiative heating and cooling of the gas.
C=======================================================================
        LOGICAL radiate
        REAL aheatmh,bheatmh2,ccoolmh2,derad,eradiate,trad,mintemp,
     &       comptmh,slowcool,meanmwt,mhboltz,fhydrogn,fhydrog2,
     &       msolunit,kpcunit

        COMMON/radcom/radiate,aheatmh,bheatmh2,ccoolmh2,eradiate,
     &                derad(nsphmax),trad,mintemp,comptmh,slowcool,
     &                meanmwt,mhboltz,fhydrogn,fhydrog2,msolunit,
     &                kpcunit

C=======================================================================
C   Definitions specific to vectorized tree construction, vectorized
C   tree walk, and vectorized tree search for nearest neighbors
C=======================================================================
        INTEGER pgroupb,subindex,bodlist,groups,templist,parent,nnear,
     &          asubp,celllist,groupbod,pnear,subpvect(nsubcell*ncells),
     &          ingroup,maxnearb,ngroups,npercell,nearbods,nnearlis,
     &          nntot,nnmin,nnmax,nnavg,isubset,nworkvec,npart
        REAL bottom,tempvect,workvect

        PARAMETER(ingroup=32,maxnearb=24*nsphmax,nworkvec=32768,
     &            npart=2*nbodsmax)

        EQUIVALENCE (subpvect(1),subp(nbods1,1))

        COMMON/nearcom/pgroupb(ncells+1),ngroups,groups(ncells),
     &                 groupbod(nbodsmax),npercell(nbods1:nbodcell),
     &                 nearbods(maxnearb),nnear(nbodsmax),
     &                 pnear(nsphmax),tempvect(nbodsmax),
     &                 nnearlis(nsphmax),bottom(1:nbodcell,ndim)
        COMMON/neighcom/nntot,nnmin,nnmax,nnavg
        COMMON/concom/celllist(ncells),parent(nbodsmax),asubp(npart),
     &                templist(nbodsmax),bodlist(nbodsmax),
     &                isubset(npart),subindex(nbodsmax)
        COMMON/workcom/workvect(nworkvec)

C=======================================================================
C   Definitions specific cosmological simulations.
C=======================================================================
        INTEGER noutcosm,maxnoutc
        LOGICAL cosmo,comove
        REAL omega0,hubble0,hubble,expanpar,redshift,exparout,ecosm,
     &       pboxsize,hboxsize,cosmofac,cosmof3,cosmohub,cosmfold

        PARAMETER(maxnoutc=100)

        COMMON/cosmcom/omega0,hubble0,hubble,expanpar,redshift,cosmo,
     &                 comove,noutcosm,exparout(maxnoutc),ecosm,
     &                 pboxsize,hboxsize,cosmofac,cosmof3,cosmohub,
     &                 cosmfold

C=======================================================================
C   Definitions specific to periodic boundary conditions.
C=======================================================================
        INTEGER nreplica,nmaxew,maxhew,nmaxew2,nmaxew2s,nmaxew2c
        LOGICAL bwrap
        REAL alphaew,hmaxew2,pi,sqrtpi,radew2,bcfx,bcfy,bcfz,bcphi,
     &       bcfxvec,bcfyvec,bcfzvec,bcphivec
     
        PARAMETER(nmaxew=32,nmaxew2=nmaxew+2,nmaxew2s=nmaxew2*nmaxew2,
     &    nmaxew2c=nmaxew2*nmaxew2s,alphaew=2.0,nreplica=3,hmaxew2=10.0)
     
        DIMENSION bcfxvec(0:nmaxew2c-1),bcfyvec(0:nmaxew2c-1),
     &    bcfzvec(0:nmaxew2c-1),bcphivec(0:nmaxew2c-1)
     
        EQUIVALENCE (bcfxvec(0),bcfx(0,0,0)),(bcfyvec(0),bcfy(0,0,0)),
     &    (bcfzvec(0),bcfz(0,0,0)),(bcphivec(0),bcphi(0,0,0))
     
        COMMON/bcparam/pi,sqrtpi,radew2,maxhew,bwrap
        COMMON/bcommon/bcfx(0:nmaxew+1,0:nmaxew+1,0:nmaxew+1),
     &            bcfy(0:nmaxew+1,0:nmaxew+1,0:nmaxew+1),
     &            bcfz(0:nmaxew+1,0:nmaxew+1,0:nmaxew+1),
     &            bcphi(0:nmaxew+1,0:nmaxew+1,0:nmaxew+1)
