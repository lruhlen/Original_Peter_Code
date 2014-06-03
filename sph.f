C***********************************************************************
C
C
                             PROGRAM treesph
C
C
C***********************************************************************
C
C
C     A code to evolve self-gravitating fluid flows using smoothed
C     particle hydrodynamics and the hierarchical tree method.  
C     The code is written in standard FORTRAN. 
C     Nearest neighbor searches are performed once per step, using a
C     tree.
C
C     In this version, vectorization of the tree walks is achieved by 
C     simultaneously processing all cells at the same level in the 
C     tree.  The gravitational force calculation and nearest neighbor 
C     search proceed for a single particle at a time, in serial order.
C
C     Smoothing is performed using a cubic spline kernel.  The kernel
C     is computed by linear interpolation from a look-up table.  The 
C     gravitational field is also smoothed according to the spline
C     kernel; however, the gravitational softening parameter is,
C     in general, not required to be equal to the SPH smoothing 
C     length.  
C     
C     In this version the SPH smoothing length is allowed to be 
C     different for each particle and is determined from the local 
C     smoothed estimate of the density.  Thus, h effectively defines
C     a sampling volume for each particle.  The smoothing lengths
C     are updated by requiring that each SPH particle interact
C     roughly with a fixed number of nearest neighbors.
C
C     Several different forms of the artificial viscosity are
C     available, including a bulk viscosity, a standard form for
C     the SPH artificial viscosity, and a version of the SPH
C     viscosity, modified by the curl of the velocity field to
C     reduce the shear component.
C
C     Periodic and quasi-periodic boundary conditions, implemented
C     according to Bouchet & Hernquist and Bouchet, Hernquist, &
C     Suto are optional.
C
C     The gravitational softening length of each particle is
C     optionally variable.  In the case of SPH particles, the
C     softening length is simply equal to the variable smoothing
C     length, whereas for collisionless particles softening lengths
C     are determined from the requirement that each softening volume
C     contain a constant number of near neighbors.  The softening
C     length for cells is always taken to be a mass-weighted mean 
C     of the softening lengths for all the particles it contains.
C
C     The computational system of units is determined by the input
C     data, with the assumption that G=1 .  Particles are not
C     required to have identical masses.
C
C
C                       Version 10: November 7, 1988
C
C
C             Lars Hernquist, Institute for Advanced Study
C             Neal Katz, University of Arizona
C
C
C=======================================================================
C
C
C     This is the top-level evolution program treesph.  Its tasks are:
C
C          1) to initialize file structures and global variables;
C          2) to input parameters and the initial system state;
C          3) to advance the state of the system for a given number
C             of timesteps;
C          4) to perform a diagnostic analysis of the system at
C             each time step (energy, angular momentum, etc.);
C          5) to periodically record the state of the system;
C          6) and to terminate the simulation and close data files.
C
C
C=======================================================================
C
C
C     Basic global variables/parameters:
C
C          acc         : acceleration components of a body.
C          atriax      : parameter defining axis ratios of triaxial
C                        mass model.
C          barmass     : mass of optional rotating, rigid bar.
C          barsize     : radius of optional rotating, rigid bar major
C                        axis.
C          boundary    : option to select boundary conditions.  Must
C                        be one of vacuum, periodic, or qperiodic.
C          btriax      : parameter defining axis ratios of triaxial
C                        mass model.
C          cellsize    : linear sizes of cells.
C          cputime     : cpu time (secs) used during the simulation.
C          cputime0    : cumulative cpu time at start of run.
C          cputime1    : cumulative cpu time at end of run.
C          ctriax      : parameter defining axis ratios of triaxial
C                        mass model.
C          dgravsum    : option to compute gravity by direct 
C                        summation (.TRUE.).
C          dtime       : the timestep.
C          dtime2      : timestep/2.
C          eps         : gravitational smoothing parameter.
C          epsvect     : values of gravitational softening length
C                        for each particle and each cell.
C          ektot       : total system kinetic energy.
C          eptot       : total system gravitational potential energy.
C          esofttot    : energy correction to virial theorem from
C                        particle softening.
C          etot        : total energy of the system.
C          fixthalo    : option to include a fixed halo contribution
C                        to the potential.
C          fixtriax    : option to include fixed triaxial mass model.
C          four        : the constant 4.
C          gamhalo     : scale-length of optional, fixed isothermal
C                        halo.
C          headline    : identification string for the run.
C          incells     : number of cells currently in use.
C          incellsg    : number of cells used to compute gravity.
C          inpteps     : option to read in values of gravitational
C                        softening parameter from particle data
C                        file -- for collisionless particles only.
C          intfour     : the integer constant 4.
C          intone      : the integer constant 1.
C          intthree    : the integer constant 3.
C          inttwo      : the integer constant 2.
C          intzero     : the integer constant 0.
C          log2        : the constant log10(2).
C          mass        : masses of bodies and cells.
C          mastriax    : mass of triaxial mass model.
C          mgas        : total mass of the gas in the system.
C          minusone    : the constant -1.
C          minustwo    : the constant -2.
C          mstar       : total mass of newly formed stars.
C          mtot        : total mass of the system.
C          nbodies     : total number of bodies.
C          nbodsmax    : maximum number of bodies.
C          ncells      : maximum number of cells.
C          ndim        : number of spatial dimensions.
C          noutbod     : frequency of system record outputs.
C          noutlog     : frequency of outputs to log file.
C          nsat        : number of particles in satellite body (used
C                        only if optional fixed halo is included).
C          nsavg       : average particle number in softening volume.
C          nsmax       : maximum particle number in softening volume.
C          nsmin       : minimum particle number in softening volume.
C          nstot       : average particle number in softening volume.
C          nsteps      : total number of timesteps in the simulation.
C          nsubcell    : number of subcells per cell.
C          nsvolume    : number of collisionless particles per
C                        softening volume; used only if variabls
C                        is .TRUE.
C          nsvtol      : fractional tolerance in number of neighbors
C                        relative to nsvolume.  A real number 
C                        typically ~ 0.05.
C          ntavg       : average length of interaction lists.
C          ntmax       : largest interaction list in current time step.
C          ntmin       : shortest interaction list in current time step.
C          nttot       : sum of interaction lists in current time step.
C          one         : the constant 1.
C          onehalf     : the constant 0.5.
C          patspeed    : pattern speed of optional rotating, rigid
C                        bar potential.
C          phi         : gravitational potential.
C          phiext      : contribution to potential from external fields.
C          pos         : coordinates of bodies, c.m. coords. of cells.
C          quad        : quadrupole moments of cells.
C          restart     : option to restart run from a SYSDUMP file.
C          rhalo       : maximum extent of optional fixed halo.
C          rigidbar    : option to include rotating, rigid bar in 
C                        force calculation.
C          rmin        : coords. of lower-left corner of system box.
C          root        : pointer to the top of the tree.
C          rsize       : length of the system box.
C          selfgrav    : option to turn off (.FALSE.) system self-
C                        gravity.
C          snheat      : total heating due to supernovae.
C          subp        : pointers to descendents of a cell.
C          tdiebar     : time over which the bar diminishes to zero
C                        strength, reckoned from the starting time
C                        plus tgrowbar plus tlivebar.
C          tgrowbar    : time over which bar grows to full strength,
C                        reckoned from the starting time.
C          three       : the constant 3.
C          tiny        : a small number used to prevent divergences.
C          tlivebar    : time over which bar is at full strength,
C                        reckoned from the starting time plus tgrowbar.
C          tnow        : current system time.
C          tol         : accuracy parameter.
C          tol2        : tol * tol.
C          tpos        : current position time.
C          ttree       : time of last update of gravitational tree.
C          two         : the constant 2.
C          usequad     : option to use (.TRUE.) quadrupole terms.
C          usesph      : option to use (.TRUE.) SPH calculation.
C          variabls    : option to use (.TRUE.) variable gravitational
C                        softening lengths for collisionless particles.
C          vel         : velocity components of a body.
C          vhalo       : asymptotic rotation velocity of optional
C                        fixed halo.
C          zero        : the constant 0.
C          zero02      : the constant 0.02.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to input/output.
C
C          ireclog                        : log file record counter.
C          uterm, upars, ulog, ubodsin,   : logical i/o unit numbers.
C            ubodsout,uboddump,utermfil,
C            ucrash,uindump
C          parsfile, logfile, ibodfile,   : character names of files.
C            obodfile,dumpfile,termfile,
C            crashfil,indumpf
C
C-----------------------------------------------------------------------
C
C   Definitions specific to individual particle time steps.
C
C          endstep     : indicator that a full step is complete.
C          etol        : tolerance parameter to select new time step.
C          inittbin    : initial time step bin of all particles.          
C          itimestp    : defines time step of each particle by 
C                        dtime/itimestp.
C          mintstep    : defines minimum allowed timestep by 
C                        dtime/mintstep.
C          npactive    : number of particles needing acceleration.
C          nsphact     : number of SPH particles needing acceleration.
C          ntvector    : minimum number of particles in each time
C                        step bin.
C          otimestp    : previous value of itimestp.
C          pactive     : indices of particles requiring acceleration.
C          stime       : fraction of large timestep remaining.
C          tsteppos    : time step with which to step positions.
C          upbin       : last time step bin which was updated.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to SPH.
C
C          acsmooth    : table of smoothed gravitational acceleration.
C          alpha, beta : coefficients in the artificial viscosity.
C          artfvisc    : option to select bulk ('bulk'), standard SPH
C                        ('sph'), or modified SPH artificial viscosity
C                        ('sphv').
C          cfrict      : value of frictional damping constant; used
C                        only if friction is .TRUE.
C          consthsm    : smoothing length if variablh is .FALSE.  A
C                        value is computed if consthsm is input zero.   
C          courant     : the courant number.
C          csound      : speed of sound of dissipative particles.
C          deldr2i     : separation of values in wsmooth.
C          dentdt      : time derivative of the entropic function a(s).
C          dentold     : time derivative of the entropic function a(s)
C                        from the previous time step.
C          dethdt      : time derivative of the thermal energy.
C          dethold     : time derivative of the thermal energy from
C                        the previous time step.
C          dwsmooth    : table of kernel derivatives.
C          entinit     : option to initialize entropic function from 
C                        a Jeans mass criterion.
C          entold      : entropic function a(s) from previous time 
C                        step.
C          entropy     : the entropic function a(s) for each particle.
C          epsfiltr    : small, dimensionless parameter used in
C                        filtering velocities of SPH particles.
C          epsgas      : gravitational smoothing parameter for gas
C                        particles.
C          epssph      : parameter in artificial viscosity.
C          ethermal    : the specific thermal energy per particle.
C          ethinit     : option to initialize thermal energy from a
C                        Jeans mass criterion.
C          ethold      : thermal energy from previous time step.
C          ethtot      : total thermal energy of the system.
C          friclucy    : option to use Lucy's method for frictional
C                        damping; otherwise a standard frictional term
C                        is used in the equations of motion.
C          friction    : option to include (.TRUE.) fricitonal damping
C                        in equations of motion for SPH particles.
C          gamma       : ratio of specific heats.
C          gamma1      : gamma - 1.
C          gamma2      : gamma - 2.
C          geogradp    : option to use geometric mean form for the
C                        pressure gradient term (.TRUE.); if .FALSE.,
C                        the arithmetic mean form is used instead.
C          hsmcurlv    : h* | curl v | for each particle.
C          hsmdivv     : h*(div v) for each particle.
C          hsmooth     : the smoothing length.
C          inpmetal    : option to input metallicity of gas and stars.
C          inptform    : option to input formation time of stars.
C          isotherm    : option to select isothermal eq. of state.
C          metals      : metallicity for gas and star particles.
C          mumaxdvh    : maximum muij in SPH viscosity or h*ABS(div v)
C                        for bulk viscosity.
C          ninterp     : number of values in look-up tables.
C          nsmooth     : hsmooth is updated so that each particle
C                        interacts with nsmooth neighbors.
C          nsmtol      : fractional tolerance in number of neighbors
C                        relative to nsmooth.  A real number 
C                        typically ~ 0.05.
C          nsph        : number of dissipational particles.
C          nstar       : number of star particles.
C          nsphmax     : maximum number of dissipational particles.
C          outputas    : option to output entropic function a(s),
C                        defined by p = a(s) * rho**gamma, to data 
C                        file rather than the temperature (default).
C          outputet    : option to output thermal energy to data
C                        file rather than the temperature (default).
C          phsmooth    : table of look-up values for grav. potential.
C          piinv       : 1.0 / 3.141592654
C          readinas    : option to read in entropic function a(s),
C                        defined by p = a(s) * rho**gamma, from
C                        data file rather than the thermal energy.
C          readinet    : option to read in thermal energy from data
C                        file rather than the temperature (default).
C                        This parameter should generally be .TRUE.
C                        if a dimensionless systems of units is used.
C          rho         : the density.
C          sphfiltr    : option to apply a spatial filter (.TRUE.)
C                        to velocities of SPH particles.
C          sphinit     : option to initialize (.TRUE.) rho and hsmooth
C                        from input particle positions and nsmooth.
C          starform    : option to include (.TRUE.) star formation.
C          symmetry    : option to select Hernquist-Katz ('hk') or
C                        Benz-Evrard ('be') symmetrized form for 
C                        smoothed estimates.  The former is defined by
C                        ~ 0.5 *(W(h1)+W(h2)), while the latter is
C                        W(0.5*(h1+h2)).
C          tempsphv    : temporary storage array for manipulating
C                        SPH data.
C          teth        : current thermal energy time.
C          tform       : formation time of stars.
C          thirty5     : the constant 35.
C          tvel        : current velocity times for SPH particles.
C          uentropy    : option to integrate the entropic function
C                        a(s) (.TRUE.) rather than the thermal energy.
C          variablg    : option to use (.TRUE.) variable gravitational
C                        softening for SPH particles, using hsmooth.
C          variablh    : option to select (.TRUE.) variable h.
C          veltpos     : velocities, time-synchronized with positions.
C          wsmooth     : table of look-up values of smoothing kernel.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to radiative heating and cooling of the gas.
C
C          aheatmh     : the constant a appearing in the heating rate
C                        equation, divided by the mass of a hydrogen
C                        atom, in c.g.s. units.
C          bheatmh2    : the constant b appearing in the heating rate
C                        equation, divided by the square of the mass of
C                        a hydrogen atom, in c.g.s. units.
C          ccoolmh2    : the constant c appearing in the cooling rate
C                        equation, divided by the square of the mass of
C                        a hydrogen atom, in c.g.s. units.
C          comptmh     : the constant appearing in the Compton cooling
C                        rate equation, divided by the mass of a
C                        hydrogen atom, in c.g.s. units.
C          ctcutoff    : cut-off factor in cooling equation.  For most
C                        applications it will lie in the range 0.3 -
C                        3.0.  (See stepeth.f.)
C          derad       : rate of change of specific thermal energy 
C                        due to radiative heating and cooling.
C          eradiate    : total, accumulated thermal energy which has
C                        been radiated from the system.
C          fhydrogn    : fraction of gas by weight in hydrogen.
C          fhydrog2    : fhydrogn ** 2.
C          kpcunit     : number of kpc per unit length.
C          meanmwt     : mean molecular weight of the gas.
C          mhboltz     : hydrogen mass divided by Boltzmann's
C                        constant, in c.g.s. units.
C          mintemp     : minimum temperature where radiative cooling
C                        is allowed.
C          msolunit    : number of solar masses per unit mass.
C          radiate     : option to include (.TRUE.) radiative effects.
C          slowcool    : maximum fractional cange in a particle's 
C                        thermal energy due to radiative and compton
C                        cooling.
C          trad        : current radiation energy time.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to vectorized tree construction, vectorized 
C   tree walk, and vectorized tree search for nearest neighbors.
C
C          asubp       : subpointers for active bodies or cells.
C          bodlist     : list of active bodies (i.e. not yet leaves).
C          bottom      : coordinates of bottom edges of cells.
C          celllist    : list of cells.
C          groupbod    : list of bodies in search groups.
C          groups      : list of grouped cells.
C          ingroup     : number of particles in groups for searching.
C          isubset     : indices of subsets of active bodies or cells.
C          maxnearb    : maximum total number of neighbors allowed.
C          nearbods    : list of neighbors for all SPH particles.
C          ngroups     : number of cells containing ingroup particles.
C          npercell    : number of particles in a cell.
C          nnavg       : average number of near neighbors.
C          nnmax       : maximum number of near neighbors.
C          nnmin       : minimum number of near neighbors.
C          nnear       : number of neighbors for SPH particles and for
C                        collisionless particles if variabls is .TRUE.
C          nnearlis    : number of neighbors in neighbor lists.
C          nntot       : total number of near neighbors.
C          nworkvec    : length of temporary work array workvect.  It
C                        should be set to MAX(9*max length of grav
C                        interaction list, max no. SPH neighbors).
C          parent      : parents of active bodies or cells.
C          pgroupb     : pointer to list of bodies in each search group.
C          pnear       : pointer to list of nearest neighbors.
C          subindex    : subindices for active bodies or cells.
C          subpvect    : vector equivalenced to subp.
C          templist    : temporary vector to swap arrays.
C          tempvect    : temporary storage vector.
C          workvect    : temporary work array.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to optional rigid black-hole potential.
C
C          amvecbh     : internal angular momentum of particle system.
C          cmposbh     : internal center of mass of particles.
C          cmvelbh     : internal center of mass velocities of 
C                        particles.
C          etabh       : strength parameter for black-hole collisions,
C                        defined by etabh = SQRT[(r_p**3/massbh)],
C                        where r_p is the pericenter distance.
C          massbh      : mass of black hole if usebh is .TRUE.
C          posbh       : relative coordinates of optional black hole.
C          tstartbh    : beginning time of black-hole encounter;
C                        usually < 0.
C          usebh       : option to include rigid black-hole potential.
C          velbh       : relative velocities of optional black hole.
C
C
C-----------------------------------------------------------------------
C
C   Definitions specific to cosmological simulations.
C
C          comove      : option to use (.TRUE.) comoving coordinates.
C          cosmfold    : value of cosmofac at previous time-step.
C          cosmo       : option to have (.TRUE.) a cosmological
C                        simulation.
C          cosmofac    : parameter set equal to the expansion parameter
C                        if cosmo and comove are true; one otherwise.
C          cosmof3     : cosmofac **3.
C          cosmohub    : cosmofac * hubble.
C          ecosm       : energy lost through cosmic expansion, if 
C                        comoving coordinates are used.
C          expanpar    : current value of the cosmological expansion
C                        parameter.
C          exparout    : values of expansion parameter at which to
C                        output system state.
C          hboxsize    : 0.5 * pboxsize.
C          hubble      : current value of the hubble constant, in
C                        appropriate units.
C          hubble0     : today's value of the hubble constant, in
C                        appropriate units.
C          maxnoutc    : maximum number of system outputs allowed.
C          noutcosm    : number of system outputs requested.
C          omega0      : for non-periodic boundaries, today's value of 
C                        the cosmological constant omega, otherwise, the
C                        value of omega at the initial time.
C          pboxsize    : box size for periodic cosmological simulations.
C          redshift    : current value of the cosmological redshift
C                        parameter.
C
C-----------------------------------------------------------------------
C
C   Definitions specific to periodic boundary conditions.
C
C          alphaew        : parameter specific to Ewald summation.
C          bcfx(nx,ny,nz) : correction force component at (0,0,0) due
C                           to a particle at (nx,ny,nz)*(1/(2*nmaxew))
C                           (direct force is not included).
C          bcfxvec        : vector equivalenced to bcfx.
C          bcfy(nx,ny,nz) : correction force y-component.
C          bcfyvec        : vector equivalenced to bcfy.
C          bcfz(nx,ny,nz) : correction force z-component.
C          bcfzvec        : vector equivalenced to bcfz.
C          bcphi(nx,ny,nz): correction potential energy.
C          bcphivec       : vector equivalenced to bcphi.
C          bwrap          : logical check for boundaries to be
C                           wrapped.
C          hmaxew2        : square of maximum radius of Ewald summation
C                           in reciprocal space.
C          maxhew         : integerized maximum radius of Ewald sum. 
C          nmaxew         : number of grid points for correction fields.
C          nmaxew2        : nmaxew + 2.
C          nmaxew2s       : nmaxew2 ** 2.
C          nmaxew2c       : nmaxew2 ** 3.
C          nreplica       : replica number in real space.  Contribution
C                           from replica particles within radius r_i,jn 
C                           < (nreplica + 0.6) * boxsize is summed.
C          pi             : the constant 3.141592654...
C          radew2         : (nreplica + 0.6) **2.
C          sqrtpi         : square root of pi.
C
C
C=======================================================================
C
C
C     Data structure used to compute gravitational field:
C
C          The body/cell tree structure is assumed to be of the
C          form discussed by Barnes and Hut.  Schematically, for
C          three dimensions (i.e. eight subcells per cell):
C
C         +-------------------------------------------------+
C  root-->| CELL:  mass, pos, quad, /, o, /, /, /, /, o, /  |
C         +----------------------------|--------------|-----+
C                                      |              |
C     +--------------------------------+              |
C     |                                               |
C     |   +----------------------------------+        |
C     +-->| BODY:  mass, pos, vel, acc, phi  |        |
C         +----------------------------------+        |
C                                                     |
C     +-----------------------------------------------+
C     |
C     |   +--------------------------------------------------+
C     +-->| CELL:  mass, pos, quad,  o, /, /, o, /, /, o, /  |
C         +--------------------------|--------|--------|-----+
C                                    |        |        |
C                                   etc.     etc.     etc.
C
C
C          The body/cell information is stored in arrays which
C          incorporate both bodies and cells.  For physical
C          quantities relevant to both bodies and cells, such as
C          mass and position, the array indices range from
C          1 --> nbodsmax + ncells.  For those quantities defined
C          only for bodies, such as velocity, the array indices
C          range from 1 --> nbodsmax.  For information relevant
C          to cells alone, such as pointers to descendants, the
C          array indices range from nbodsmax + 1 --> nbodsmax +
C          ncells.  With this convention, pointers can refer to
C          either bodies or cells without conflict.
C
C          The identification of a given unit as a body or a cell
C          is determined by the pointer to the body/cell.  For a
C          body, p is less than or equal to nbodsmax, while for a
C          cell, p > nbodsmax.
C
C          For applications using the SPH formalism, the dissipative
C          particles are assumed to occupy the first nsph slots.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================

C   Initialize state of the system.
C   -------------------------------
         write(*,*) 'start'
        CALL initsys
C            -------

C   Advance system state for a given number of steps.
C   -------------------------------------------------

        DO 100 n=1,nsteps

           CALL stepsys(n)
C               -------

 100    CONTINUE

C   Terminate the simulation.
C   -------------------------
        CALL endrun
C            ------

        STOP
        END
C***********************************************************************
C
C
                          SUBROUTINE initsys
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the state of the system.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        REAL second

C=======================================================================

C   Begin timing.
C   -------------
        cputime0=SECOND()

C-----------------------------------------------------------------------
C   Open data files, read input parameters and initial system state,
C   and initialize system parameters.
C-----------------------------------------------------------------------
        CALL startout
C            --------
        CALL inparams
C            --------

        IF(.NOT.restart) THEN

           CALL inbods
C               ------
        ELSE

           CALL indump
C               ------
        ENDIF

        CALL checkinp
C            --------
        CALL initpars
C            --------
        IF(boundary.EQ.'periodic') CALL initbc
C                                       ------

        IF(.NOT.restart) CALL initstep
C                             --------

        IF(.NOT.restart) THEN

C   Zero out potential and acceleration.
C   ------------------------------------
           CALL zeroacc
C               -------
           CALL zeropot
C               -------

C   Initialize hydrodynamic properties.
C   -----------------------------------
           IF(usesph.AND.nsph.GT.0) CALL initsph
C                                        -------

C   Initialize gravitational softening lengths.
C   -------------------------------------------
           IF(variabls.OR.eps.EQ.0.0) CALL maketeps
C                                          --------
           IF(.NOT.inpteps) CALL initeps
C                                -------
           IF(variabls) CALL neighcol('correct')
C                            --------

C   Compute gravitational potential and acceleration.
C   -------------------------------------------------
           CALL gravity('both')
C               -------

        ENDIF

C   Output system state.
C   --------------------
        IF(.NOT.restart) THEN
           CALL outstate(0)
C               --------
        ELSE
           CALL outhead
C               -------
        ENDIF

        IF(.NOT.restart) CALL initpos
C                             -------

        IF(restart) THEN

           CALL setbox('sph ')
C               ------
           CALL loadtree('sph ')
C               --------
           CALL cellnumb
C               --------
           CALL celledge
C               --------
           CALL groupcel
C               --------
           CALL neighbor('predict')
C               --------
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE startout
C
C
C***********************************************************************
C
C
C     Subroutine to open disk files for subsequent input/output.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Open log file.
C   --------------
        OPEN(UNIT=ulog,FILE=logfile,STATUS='NEW',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

        ireclog=0

        WRITE(ulog,10,REC=1)
 10     FORMAT(' Start of logfile output ')

        CLOSE(UNIT=ulog)
 
C   Create terminal emulation file.
C   -------------------------------
        OPEN(UNIT=utermfil,FILE=termfile,STATUS='UNKNOWN')
        WRITE(utermfil,*) ' Start of output '
        CLOSE(UNIT=utermfil)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE inparams
C
C
C***********************************************************************
C
C
C     Subroutine to read in parameters.
C
C     Input parameters:
C
C        headline  : identification string for the run.
C        nsteps    : number of timesteps.
C        noutbod   : output system state once every nsteps/noutbod 
C                    steps.
C        noutlog   : output logfile data once every nsteps/noutlog
C                    steps.
C        dtime     : the timestep.
C        tol       : error tolerance; 0.0 => exact (PP) calculation.
C        eps       : potential softening parameter.
C        inpteps   : option to read gravitational softening lengths
C                    of collisionless particles from input data.
C        variabls  : option to use variable gravitational softening
C                    lengths for collisionless particles.
C        nsvolume  : number of collisionless particles per softening
C                    volume.
C        nsvtol    : fractional tolerance in number of neighbors
C                    relative to nsvolume; typically ~ 0.05.
C        usequad   : option to include (.TRUE.) quadrupole terms.
C        dgravsum  : option to use direct summation gravity (.TRUE.).
C        boundary  : parameter to select boundary condition -- one
C                    of vacuum, periodic, or qperiodic.
C        restart   : option to restart run from a SYSDUMP file.
C-----------------------------------------------------------------------
C        inittbin  : initial time step bin for all particles.
C        mintstep  : parameter defining minimum allowed time step.
C        etol      : energy tolerance to select timestep.
C        ntvector  : minumum number of particles per time step bin.
C-----------------------------------------------------------------------
C        usesph    : option to use (.TRUE.) an SPH calculation.
C        sphinit   : option to initialize rho and hsmooth.
C        uentropy  : option to integrate entropic function a(s)
C                    (.TRUE.) rather than the thermal energy.
C        isotherm  : option to use isothermal eq. of state.
C        artfvisc  : artificial viscosity option parameter.
C        epsgas    : potential softening parameter for the gas
C                    particles.
C        gamma     : ratio of specific heats.
C        alpha     : coefficient in artificial viscosity.
C        beta      : coefficient in artificial viscosity.
C        epssph    : parameter in artificial viscosity.
C        courant   : courant number.
C        variablh  : option to have variable smoothing lengths.
C        variablg  : option to have variable SPH softening lengths.
C        consthsm  : smoothing length if variablh is .FALSE.  A value 
C                    of consthsm is computed if it is zero initially.
C        nsmooth   : number of neighbors to interact with.
C        nsmtol    : fractional tolerance in number of neighbors
C                    relative to nsmooth; typically ~ 0.05.
C        symmetry  : option to select symmetrized form of smoothed
C                    estimates.
C        geogradp  : option to use geometric mean form for pressure
C                    gradient term; if .FALSE., arithmetic mean used.
C        sphfiltr  : option to apply a filter to velocities of SPH
C                    particles.
C        epsfiltr  : small parameter used in filtering SPH velocities.
C        readinas  : option to input entropic function, a(s), rather
C                    than temperature.
C        outputas  : option to output entropic function, a(s), rather
C                    than temperature.
C        readinet  : option to input thermal energy rather than
C                    temperature.
C        outputet  : option to output thermal energy rather than
C                    temperature.
C        inpmetal  : option to input metallicity for gas and star
C                    particles.
C        inptform  : option to input formation time of stars.
C-----------------------------------------------------------------------
C        radiate   : option to include (.TRUE.) radiative effects.
C        aheatmh   : heating constant a / hydrogen mass.
C        bheatmh2  : heating constant b / hydrogen mass **2.
C        ccoolmh2  : cooling constant c / hydrogen mass **2.
C        comptmh   : Compton cooling constant / hydrogen mass.
C        meanmwt   : mean molecular weight.
C        fhydrogn  : fraction of gas by weight in hydrogen.
C        mhboltz   : hydrogen mass / Boltzmann constant.
C        kpcunit   : number of kpc per unit length.
C        msolunit  : number of solar masses per unit length.
C        mintemp   : minimum temperature for radiative cooling.
C        slowcool  : maximum fractional change in thermal energy due
C                    to cooling.
C        ctcutoff  : coefficient in cooling cut-off expression.
C        ethinit   : option to initialize (.TRUE.) thermal energy
C                    from a Jeans mass criterion. 
C        entinit   : option to initialize (.TRUE.) entropic function
C                    a(s) from a Jeans mass criterion. 
C-----------------------------------------------------------------------
C        fixthalo  : option to have a include a fixed halo potential.
C        vhalo     : asymptotic rotation velocity of fixed halo.
C        gamhalo   : scale-length of fixed halo.
C        rhalo     : maximum extent of fixed halo.
C        nsat      : number of particles in satellite body.
C-----------------------------------------------------------------------
C        fixtriax  : option to include fixed triaxial mass model.
C        atriax    : parameter defining axis ratios of triaxial model.
C        btriax    : parameter defining axis ratios of triaxial model.
C        ctriax    : parameter defining axis ratios of triaxial model.
C        mastriax  : mass of triaxial mass model.
C-----------------------------------------------------------------------
C        rigidbar  : option to include rotating, rigid bar potential.
C        patspeed  : pattern speed of rotating, rigid bar.
C        barmass   : mass of rotating, rigid bar.
C        barsize   : radius of rotating, rigid bar major axis.
C        tgrowbar  : time over which bar grows to full strength.
C        tlivebar  : time over which bar is at full strength.
C        tdiebar   : time over which bar diminishes to zero strength.
C-----------------------------------------------------------------------
C        usebh     : option to include rigid black-hole potential.
C        massbh    : mass of optional black hole.
C        etabh     : strength parameter for black-hole collisions.
C        tstartbh  : beginning time for black-hole encounter.
C-----------------------------------------------------------------------
C        selfgrav  : option to turn off (.FALSE.) system self-gravity.
C        starform  : option to include star formation.
C        friction  : option to include frictional damping for SPH
C                    particles.
C        friclucy  : option to use Lucy's damping method; otherwise a
C                    standard frictional term is included.
C        cfrict    : value of frictional damping constant.  Used only
C                    if friction is .TRUE.
C-----------------------------------------------------------------------
C        omega0    : current value of cosmological omega for vacuum
C                    b.c., initial value of omega for periodic b.c.
C        hubble0   : current value of the hubble constant.
C        cosmo     : option to have a cosmological simulation.
C        comove    : option to have comoving coordinates.
C        ecosm     : initial value of energy lost through cosmic
C                    expansion (should be zero except for restarts).
C        pboxsize  : length of periodic box.
C        noutcosm  : number of system outputs for cosmological 
C                    simulations.
C        exparout  : expansion parameters at which to output state
C                    of the system.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER *1 pcomment
        INTEGER i

C=======================================================================
 
        OPEN(UNIT=upars,FILE=parsfile,STATUS='OLD')

C   Read parameters, close the file.
C   --------------------------------

        READ(upars,'(a)') pcomment

        READ(upars,'(a)') headline
        READ(upars,*) nsteps
        READ(upars,*) noutbod
        READ(upars,*) noutlog
        READ(upars,*) dtime
        READ(upars,*) tol
        READ(upars,*) eps
        READ(upars,*) inpteps
        READ(upars,*) variabls
        READ(upars,*) nsvolume
        READ(upars,*) nsvtol
        READ(upars,*) usequad
        READ(upars,*) dgravsum
        READ(upars,'(a)') boundary
        READ(upars,*) restart

        READ(upars,'(a)') pcomment

        READ(upars,*) inittbin
        READ(upars,*) mintstep
        READ(upars,*) etol
        READ(upars,*) ntvector

        READ(upars,'(a)') pcomment

        READ(upars,*) usesph
        READ(upars,*) sphinit
        READ(upars,*) uentropy
        READ(upars,*) isotherm
        READ(upars,'(a)') artfvisc
        READ(upars,*) epsgas
        READ(upars,*) gamma
        READ(upars,*) alpha
        READ(upars,*) beta
        READ(upars,*) epssph
        READ(upars,*) courant
        READ(upars,*) variablh
        READ(upars,*) variablg
        READ(upars,*) consthsm
        READ(upars,*) nsmooth
        READ(upars,*) nsmtol
        READ(upars,'(a)') symmetry
        READ(upars,*) geogradp
        READ(upars,*) sphfiltr
        READ(upars,*) epsfiltr
        READ(upars,*) readinas
        READ(upars,*) outputas
        READ(upars,*) readinet
        READ(upars,*) outputet
        READ(upars,*) inpmetal
        READ(upars,*) inptform

        READ(upars,'(a)') pcomment

        READ(upars,*) radiate
        READ(upars,*) aheatmh
        READ(upars,*) bheatmh2
        READ(upars,*) ccoolmh2
        READ(upars,*) comptmh
        READ(upars,*) meanmwt
        READ(upars,*) fhydrogn
        READ(upars,*) mhboltz
        READ(upars,*) kpcunit
        READ(upars,*) msolunit
        READ(upars,*) mintemp
        READ(upars,*) slowcool
        READ(upars,*) ctcutoff
        READ(upars,*) ethinit
        READ(upars,*) entinit

        aheatmh=aheatmh*3.46e14*kpcunit**2.5/msolunit**1.5
        bheatmh2=bheatmh2*2.34e-17/SQRT(kpcunit*msolunit)
        ccoolmh2=ccoolmh2*2.34e-17/SQRT(kpcunit*msolunit)
        mhboltz=mhboltz*4.30e4*msolunit/kpcunit
        comptmh=comptmh*3.46e14*kpcunit**2.5/msolunit**1.5

        READ(upars,'(a)') pcomment

        READ(upars,*) fixthalo
        READ(upars,*) vhalo
        READ(upars,*) gamhalo
        READ(upars,*) rhalo
        READ(upars,*) nsat

        READ(upars,'(a)') pcomment

        READ(upars,*) fixtriax
        READ(upars,*) atriax
        READ(upars,*) btriax
        READ(upars,*) ctriax
        READ(upars,*) mastriax

        READ(upars,'(a)') pcomment

        READ(upars,*) rigidbar
        READ(upars,*) patspeed
        READ(upars,*) barmass
        READ(upars,*) barsize
        READ(upars,*) tgrowbar
        READ(upars,*) tlivebar
        READ(upars,*) tdiebar

        READ(upars,'(a)') pcomment

        READ(upars,*) usebh
        READ(upars,*) massbh
        READ(upars,*) etabh
        READ(upars,*) tstartbh

        READ(upars,'(a)') pcomment

        READ(upars,*) selfgrav
        READ(upars,*) starform
        READ(upars,*) friction
        READ(upars,*) friclucy
        READ(upars,*) cfrict

        READ(upars,'(a)') pcomment

        READ(upars,*) omega0
        READ(upars,*) hubble0
        READ(upars,*) cosmo
        READ(upars,*) comove
        READ(upars,*) ecosm
        READ(upars,*) pboxsize
        READ(upars,*) noutcosm

        IF(noutcosm.GE.maxnoutc) CALL terror(' i/o error in inparams')

        DO 10 i=1,noutcosm
           READ(upars,*) exparout(i)
 10     CONTINUE

        DO 20 i=noutcosm+1,maxnoutc
           exparout(i)= -1.
 20     CONTINUE

        CLOSE(UNIT=upars)
 
        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE inbods
C
C
C***********************************************************************
C
C
C     Subroutine to read in the data associated with the bodies.  The
C     records are assumed to be of the form: nbodies, ndim, time,
C     mass(1)...mass(n), x(1)...x(n), y(1)...y(n), ..., vx(1)...vx(n),
C     vy(1)...vy(n), ..., followed by the SPH data, if needed.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,ndimi

C=======================================================================
 
        OPEN(UNIT=ubodsin,FILE=ibodfile,STATUS='OLD')

C   Read in body data.
C   ------------------

        READ(ubodsin,*) nbodies,nsph,nstar
        READ(ubodsin,*) ndimi
        READ(ubodsin,*) tnow

        IF(nbodies.GT.nbodsmax.OR.ndimi.NE.ndim)
     &     CALL terror(' error in inbods--inconsistent inputs ')
C               ------

        DO 10 p=1,nbodies
           READ(ubodsin,*) mass(p)
 10     CONTINUE

        DO 20 p=1,nbodies
           READ(ubodsin,*) pos(p,1)
 20     CONTINUE
        DO  p=1,nbodies
           READ(ubodsin,*) pos(p,2)
        end do
        DO p=1,nbodies
           READ(ubodsin,*) pos(p,3)
        end do

        DO 30 p=1,nbodies
           READ(ubodsin,*) vel(p,1)
 30     CONTINUE
        DO  p=1,nbodies
           READ(ubodsin,*) vel(p,2)
        end do
        DO  p=1,nbodies
           READ(ubodsin,*) vel(p,3)
        end do

        IF(inpteps) THEN
           DO 40 p=nsph+1,nbodies
              READ(ubodsin,*) epsvect(p)
 40        CONTINUE
        ENDIF

        IF(.NOT.sphinit) THEN
           DO 50 p=1,nsph
              READ(ubodsin,*) rho(p)
 50        CONTINUE
        ENDIF

        IF((.NOT.ethinit).AND.(.NOT.entinit)) THEN
           DO 60 p=1,nsph
              READ(ubodsin,*) tempvect(p)
 60        CONTINUE
        ENDIF

        IF(.NOT.sphinit) THEN
           DO 70 p=1,nsph
              READ(ubodsin,*) hsmooth(p)
 70        CONTINUE
        ENDIF

        IF(inpmetal) THEN
           DO 80 p=1,nsph
              READ(ubodsin,*) metals(p)
 80        CONTINUE

           DO 90 p=nbodies-nstar+1,nbodies
              READ(ubodsin,*) metals(p)
 90        CONTINUE
        ENDIF

        IF(inptform) THEN
           DO 100 p=nbodies-nstar+1,nbodies
              READ(ubodsin,*) tform(p)
 100       CONTINUE
        ENDIF

        IF((.NOT.ethinit).AND.(.NOT.entinit)) THEN
           IF(readinas.OR.readinet) THEN
              IF(.NOT.uentropy) THEN
                 DO 200 p=1,nsph
                    ethermal(p)=tempvect(p)
 200             CONTINUE
              ELSE
                 DO 210 p=1,nsph
                    entropy(p)=tempvect(p)
 210             CONTINUE
              ENDIF
           ELSE
              IF(.NOT.uentropy) THEN
                 IF(.NOT.isotherm) THEN
                    DO 220 p=1,nsph
                       ethermal(p)=tempvect(p)/(meanmwt*
     &                             mhboltz*(gamma-1.))
 220                CONTINUE
                 ELSE
                    DO 225 p=1,nsph
                       ethermal(p)=tempvect(p)*gamma/
     &                             (meanmwt*mhboltz)
 225                CONTINUE
                 ENDIF
              ELSE
                 DO 230 p=1,nsph
                    entropy(p)=tempvect(p)
 230             CONTINUE
              ENDIF
           ENDIF
        ENDIF

        CLOSE(ubodsin)
 
        RETURN
        END
C***********************************************************************
C
C
                            SUBROUTINE indump
C
C
C***********************************************************************
C
C
C     Subroutine to read in system dump from an ascii data file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p

C=======================================================================
 
        OPEN(UNIT=uindump,FILE=indumpf,STATUS='OLD')

C   Output system state.
C   --------------------

        READ(uindump,*) nbodies,nsph,nstar

        READ(uindump,*) tnow,tpos
        READ(uindump,*) dtime
 
        READ(uindump,*) mtot
        READ(uindump,*) ektot
        READ(uindump,*) eptot
        READ(uindump,*) snheat
        READ(uindump,*) eradiate
        READ(uindump,*) esofttot
        READ(uindump,*) enerror

        READ(uindump,*) teth
        READ(uindump,*) upbin
        READ(uindump,*) stime
        READ(uindump,*) tsteppos
        READ(uindump,*) endstep
        READ(uindump,*) trad

        DO 70 p=1,nbodies
           READ(uindump,*) mass(p)
 70     CONTINUE

        DO 80 p=1,nbodies
           READ(uindump,*) pos(p,1),pos(p,2),pos(p,3)
 80     CONTINUE

        DO 90 p=1,nbodies
           READ(uindump,*) vel(p,1),vel(p,2),vel(p,3)
 90     CONTINUE

        DO 105 p=1,nbodies
           READ(uindump,*) epsvect(p)
 105    CONTINUE

        DO 100 p=1,nsph
           READ(uindump,*) rho(p)
 100    CONTINUE

        DO 110 p=1,nsph
           READ(uindump,*) ethermal(p)
 110    CONTINUE

        DO 120 p=1,nsph
           READ(uindump,*) hsmooth(p)
 120    CONTINUE

        DO 130 p=1,nsph
           READ(uindump,*) csound(p)
 130    CONTINUE

        DO 140 p=1,nsph
           READ(uindump,*) ethold(p)
 140    CONTINUE

        DO 150 p=1,nsph
           READ(uindump,*) dethdt(p)
 150    CONTINUE

        DO 160 p=1,nsph
           READ(uindump,*) dethold(p)
 160    CONTINUE

        DO 170 p=1,nsph
           READ(uindump,*) tvel(p)
 170    CONTINUE

        DO 180 p=1,nsph
           READ(uindump,*) mumaxdvh(p)
 180    CONTINUE

        DO 190 p=1,nsph
           READ(uindump,*) hsmdivv(p)
 190    CONTINUE

        DO 200 p=1,nsph
           READ(uindump,*) derad(p)
 200    CONTINUE

        DO 204 p=1,nsph
           READ(uindump,*) hsmcurlv(p)
 204    CONTINUE

        DO 210 p=1,nbodies
           READ(uindump,*) nnear(p)
 210    CONTINUE

        DO 220 p=1,nbodies
           READ(uindump,*) phi(p),phiext(p)
 220    CONTINUE

        DO 230 p=1,nbodies
           READ(uindump,*) acc(p,1),acc(p,2),acc(p,3)
 230    CONTINUE

        DO 240 p=1,nbodies
           READ(uindump,*) itimestp(p)
 240    CONTINUE

        DO 250 p=1,nbodies
           READ(uindump,*) otimestp(p)
 250    CONTINUE

        DO 260 p=1,nsph
           READ(uindump,*) metals(p)
 260    CONTINUE

        DO 270 p=nbodies-nstar+1,nbodies
           READ(uindump,*) metals(p)
 270    CONTINUE

        DO 280 p=nbodies-nstar+1,nbodies
           READ(uindump,*) tform(p)
 280    CONTINUE

        CLOSE(UNIT=uindump)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE checkinp
C
C
C***********************************************************************
C
C
C     Subroutine to check consistency of input parameters and data,
C     output warnings to the terminal and/or log file, and terminate
C     the simulation if necessary.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

        IF(nsteps.LT.0.OR.nsteps.GT.10000000)
     &     CALL terror(' input error for parameter nsteps ')
C               ------

        IF(noutbod.LT.0)
     &     CALL terror(' input error for parameter noutbod ')
C               ------

        IF(noutlog.LT.0)
     &     CALL terror(' input error for parameter noutlog ')
C               ------

        IF(dtime.LE.0.0.OR.dtime.GT.1.e20)
     &     CALL terror(' input error for parameter dtime ')
C               ------

        IF((.NOT.dgravsum).AND.(tol.LT.0.0.OR.tol.GT.1.5))
     &     CALL terror(' input error for parameter tol ')
C               ------

        IF(eps.LT.0.0.OR.eps.GT.1.e20)
     &     CALL terror(' input error for parameter eps ')
C               ------

        IF(variabls.AND.(nsvolume.LE.0.OR.nsvolume.GT.nbodsmax))
     &     CALL terror(' input error for parameter nsvolume ')
C               ------

        IF(boundary.NE.'vacuum'.AND.boundary.NE.'periodic'.AND.
     &     boundary.NE.'qperiodic')
     &     CALL terror(' input error for parameter boundary ')
C               ------

        IF(inittbin.LT.1.OR.inittbin.GT.262144.OR.inittbin.GT.mintstep)
     &     CALL terror(' input error for parameter inittbin ')
C               ------

        IF(mintstep.LT.1.OR.mintstep.GT.262144)
     &     CALL terror(' input error for parameter mintstep ')
C               ------

        IF(etol.LE.0.0.OR.etol.GT.1.e20)
     &     CALL terror(' input error for parameter etol ')
C               ------

        IF(ntvector.LE.0.OR.ntvector.GT.262144)
     &     CALL terror(' input error for parameter ntvector ')
C               ------

        IF(usesph) THEN

           IF(artfvisc.NE.'bulk'.AND.artfvisc.NE.'sphv'.AND.
     &        artfvisc.NE.'sph ')
     &        CALL terror(' input error for parameter artfvisc ')
C                  ------

           IF(epsgas.LE.0.0.OR.epsgas.GT.1.e20)
     &        CALL terror(' input error for parameter epsgas ')
C               ------

           IF(gamma.LE.0.0.OR.gamma.GT.1.e20)
     &        CALL terror(' input error for parameter gamma ')
C                  ------

           IF(alpha.LT.0.0.OR.alpha.GT.1.e20)
     &        CALL terror(' input error for parameter alpha ')
C                  ------

           IF(beta.LT.0.0.OR.beta.GT.1.e20)
     &        CALL terror(' input error for parameter beta ')
C                  ------

           IF(epssph.LE.0.0.OR.epssph.GT.1.e20)
     &        CALL terror(' input error for parameter epssph ')
C                  ------

           IF(courant.LE.0.0.OR.courant.GT.1.e20)
     &        CALL terror(' input error for parameter courant ')
C                  ------

           IF(consthsm.LT.0.0.OR.consthsm.GT.1.e20.AND.(.NOT.variablh))
     &        CALL terror(' input error for parameter consthsm ')
C                  ------

           IF(nsmooth.LE.0.OR.nsmooth.GT.nbodsmax.AND.variablh)
     &        CALL terror(' input error for parameter nsmooth ')
C                  ------

           IF(symmetry.NE.'hk'.AND.symmetry.NE.'be')
     &        CALL terror(' input error for parameter symmetry ')
C                  ------

           IF(aheatmh.LT.0.0.OR.aheatmh.GT.1.e20.AND.radiate)
     &        CALL terror(' input error for parameter aheatmh ')
C                  ------

           IF(bheatmh2.LT.0.0.OR.bheatmh2.GT.1.e20.AND.radiate)
     &        CALL terror(' input error for parameter bheatmh2 ')
C                  ------
 
           IF(ccoolmh2.LT.0.0.OR.ccoolmh2.GT.1.e20.AND.radiate)
     &        CALL terror(' input error for parameter ccoolmh2 ')
C                  ------

           IF(comptmh.LT.0.0.OR.comptmh.GT.1.e20.AND.radiate)
     &        CALL terror(' input error for parameter comptmh ')
C                  ------

           IF(meanmwt.LE.0.0.OR.meanmwt.GT.1.e20)
     &        CALL terror(' input error for parameter meanmwt ')
C                  ------

           IF(fhydrogn.LT.0.0.OR.fhydrogn.GT.1.0)
     &        CALL terror(' input error for parameter fhydrogn ')
C                  ------

           IF(mhboltz.LE.0.0.OR.mhboltz.GT.1.e20)
     &        CALL terror(' input error for parameter mhboltz ')
C                  ------

           IF(mintemp.LT.0.0.OR.mintemp.GT.1.e20)
     &        CALL terror(' input error for parameter mintemp ')
C                  ------

           IF(slowcool.LT.0.0.OR.slowcool.GT.1.e20.AND.radiate)
     &        CALL terror(' input error for parameter slowcool ')
C                  ------

           IF(ctcutoff.LT.0.0.OR.ctcutoff.GT.1.e20.AND.radiate)
     &        CALL terror(' input error for parameter ctcutoff ')
C                  ------

        ENDIF

        IF(fixthalo) THEN

           IF(vhalo.LT.0.0.OR.vhalo.GT.1.e20)
     &        CALL terror(' input error for parameter vhalo ')
C                  ------

           IF(gamhalo.LT.0.0.OR.gamhalo.GT.1.e20)
     &        CALL terror(' input error for parameter gamhalo ')
C                  ------

           IF(rhalo.LT.0.0.OR.rhalo.GT.1.e20)
     &        CALL terror(' input error for parameter rhalo ')
C                  ------

           IF(nsat.LT.0.OR.nsat.GT.nbodsmax)
     &        CALL terror(' input error for parameter nsat ')
C                  ------

        ENDIF

        IF(fixtriax) THEN

           IF(atriax.LT.0.0.OR.atriax.GT.1.e20)
     &        CALL terror(' input error for parameter atriax ')
C                  ------

           IF(btriax.LT.0.0.OR.btriax.GT.1.e20)
     &        CALL terror(' input error for parameter btriax ')
C                  ------

           IF(ctriax.LT.0.0.OR.ctriax.GT.1.e20)
     &        CALL terror(' input error for parameter ctriax ')
C                  ------

           IF(mastriax.LT.0.OR.mastriax.GT.1.e20)
     &        CALL terror(' input error for parameter mastriax ')
C                  ------

           IF(atriax.LT.btriax.OR.atriax.LT.ctriax)
     &        CALL terror(' input error for parameter atriax ')
C                  ------

           IF(btriax.GT.atriax.OR.btriax.LT.ctriax)
     &        CALL terror(' input error for parameter btriax ')
C                  ------

           IF(ctriax.GT.atriax.OR.ctriax.GT.btriax)
     &        CALL terror(' input error for parameter ctriax ')
C                  ------

        ENDIF

        IF(usebh) THEN

           IF(massbh.LT.0.0.OR.massbh.GT.1.e20)
     &        CALL terror(' input error for parameter massbh ')
C               ------

           IF(etabh.LT.0.0.OR.etabh.GT.1.e20)
     &        CALL terror(' input error for parameter etabh ')
C               ------

           IF(tstartbh.LT.-1.e20.OR.tstartbh.GT.1.e20)
     &        CALL terror(' input error for parameter tstartbh ')
C               ------
        ENDIF

        IF(cfrict.LT.0.0.OR.cfrict.GT.1.e20.AND.friction)
     &     CALL terror(' input error for parameter cfrict ')
C               ------

        IF(cosmo) THEN

           IF(omega0.LT.0.0.OR.omega0.GT.1.e20)
     &        CALL terror(' input error for parameter omega0 ')
C                  ------

           IF(hubble0.LT.0.0.OR.hubble0.GT.1.e20)
     &        CALL terror(' input error for parameter hubble0 ')
C                  ------

           IF(ecosm.LE.-1.e20.OR.ecosm.GT.1.e20.AND.comove)
     &        CALL terror(' input error for parameter ecosm ')
C                  ------

           IF(noutcosm.LT.0.OR.noutcosm.GT.262144.AND.comove)
     &        CALL terror(' input error for parameter noutcosm ')
C                  ------

        ENDIF

        IF(friction.AND.(alpha.GT.0.0.OR.beta.GT.0.0)) THEN
           CALL outterm(' WARNING -- friction and viscosity used ',0)
C               -------
           CALL outterm(' alpha and beta set to zero ',0)
C               -------
           alpha=0.0
           beta=0.0
        ENDIF

        IF(friction.AND.radiate) THEN
           CALL outterm(' WARNING -- friction and radiate used ',0)
C               -------
           CALL outterm(' radiate set to .FALSE. ',0)
C               -------
           radiate=.FALSE.
        ENDIF

        IF(usesph.AND.nworkvec.LT.nsph)
     &     CALL terror(' array bound error in checkinp ')
C               ------

        IF((boundary.EQ.'periodic'.OR.boundary.EQ.'qperiodic').AND.
     &      pboxsize.LE.0.0.OR.pboxsize.GT.1.e20)
     &     CALL terror(' input error for parameter pboxsize ')
C               ------

        IF((.NOT.cosmo).AND.comove) 
     &     CALL terror(' comove allowed only if cosmo is true ')
C               ------

        IF(usequad.AND.boundary.NE.'vacuum')
     &     CALL terror(' quad. option allowed only with vaccum b.c. ')
C               ------

        IF(friction.AND.comove)
     &     CALL terror(' friction and comove not allowed ')
C               ------

        IF(isotherm.AND.ABS(gamma-1.0).GT.1.e-6)
     &     CALL terror(' gamma must be 1.0 for isothermal gas')
C               ------

        IF(uentropy.AND.isotherm)
     &     CALL terror(' uentropy must be .FALSE. for isothermal gas')
C               ------

        IF(ethinit.AND.isotherm)
     &     CALL terror(' ethinit must be .FALSE. for isothermal gas')
C               ------

        IF(entinit.AND.isotherm)
     &     CALL terror(' entinit must be .FALSE. for isothermal gas')
C               ------

        IF(readinet.AND.isotherm)
     &     CALL terror(' readinet must be .FALSE. for isothermal gas')
C               ------

        IF(readinas.AND.isotherm)
     &     CALL terror(' readinas must be .FALSE. for isothermal gas')
C               ------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpars
C
C
C***********************************************************************
C
C
C     Subroutine to initialize system parameters that depend on
C     either the input data or defined PARAMETERS.  The local
C     variable p is a pointer to the bodies/cells.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i
        REAL deldr2,xw,xw2,deldrg,xw3,xw4

C=======================================================================

C   Initialize misc. useful numbers.
C   --------------------------------
        minusone = -1.
        minustwo = -2.
        tiny=1.e-20
        zero=0.
        zero02=0.02
        intzero=0
        intone=1
        inttwo=2
        intthree=3
        intfour=4
        one=1.
        two=2.
        three=3.
        four=4.
        log2=LOG(two)
        onehalf=0.5

        bwrap = (boundary.EQ.'periodic').OR.(boundary.EQ.'qperiodic')

        IF(usebh.AND..NOT.restart) THEN
           tnow=tstartbh
        ENDIF

C   Initialize position and velocity times, 1/2 timestep.
C   -----------------------------------------------------
        IF(.NOT.restart) THEN
           tpos=tnow
           teth=tnow
           trad=tnow
        ENDIF

        ttree=tnow-1.
        dtime2=.5*dtime
        tol2=tol*tol

C   Initialize size parameter for bodies.
C   -------------------------------------

        DO 5 p=1,nbodies
           cellsize(p)=0.
 5      CONTINUE

        rsize=0.
        rmin(1)=0.
        rmin(2)=0.
        rmin(3)=0.

C=======================================================================
C   Initialize SPH parameters.
C=======================================================================
        gamma1 = gamma - 1.
        gamma2 = gamma - 2.
        piinv= 1.0 / (4.*ATAN(one))

        IF(.NOT.restart) THEN
           eradiate = 0.0
           snheat = 0.0
           enerror = 0.0
        ENDIF

        fhydrog2=fhydrogn**2
        thirty5=35.

        IF(.NOT.usesph) nsph=0

        IF(.NOT.restart) THEN

           DO 7 p=1,nbodies
              nnear(p)=0
 7         CONTINUE

           DO 10 p=1,nsph
              tvel(p)=tnow
              mumaxdvh(p)=0.0
              csound(p)=0.0
 10        CONTINUE

           IF(usesph.AND.sphinit) THEN
              DO 15 p=1,nsph
                 hsmooth(p)=0.0
 15           CONTINUE
           ENDIF

        ENDIF

C   Initialize cosmological parameters.
C   -----------------------------------
        IF(.NOT.cosmo) THEN
           ecosm=0.0
           comptmh=0.0
           redshift=0.0
           expanpar=1.0
        ENDIF

        cosmofac=1.0
        cosmfold=1.0
        cosmof3=1.0
        cosmohub=0.0

        hboxsize=0.5*pboxsize
        pi=ACOS(-one)
        sqrtpi=SQRT(pi)

        IF(boundary.EQ.'periodic') THEN
           radew2=(REAL(nreplica)+0.6)**2
           maxhew=INT(SQRT(hmaxew2))
        ENDIF

C-----------------------------------------------------------------------
C   Initialize variables and arrays for smoothing kernel interpolation.
C   Interpolation performed in distance **2, using a spline kernel.
C-----------------------------------------------------------------------
        deldr2=4./ninterp
        deldr2i=1./deldr2

        DO 20 i=0,1+ninterp
           xw2=i*deldr2
           xw=SQRT(xw2)
           IF(xw2.LE.one) THEN
              wsmooth(i)=one-1.5*xw2+0.75*xw*xw2
              dwsmooth(i)=-3.+2.25*xw
           ELSE
              wsmooth(i)=0.25*(2.-xw)**3
              dwsmooth(i)=-0.75*(2.-xw)**2/(xw+tiny)
           ENDIF
           IF(xw2.GE.four) THEN
              wsmooth(i)=zero
              dwsmooth(i)=zero
           ENDIF
 20     CONTINUE

C-----------------------------------------------------------------------
C   Initialize variables and arrays for gravitational field smoothing 
C   interpolation.  Interpolation performed in distance.
C-----------------------------------------------------------------------
        deldrg=2./ninterp

        phsmooth(0)=7.*SQRT(tiny)/5.
        acsmooth(0)=4.*SQRT(tiny)*tiny/3.

        DO 30 i=1,1+ninterp
           xw=i*deldrg
           xw2=xw*xw
           xw3=xw2*xw
           xw4=xw2*xw2
           IF(xw.LE.one) THEN
              phsmooth(i)=-2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
              acsmooth(i)=xw3*(4./3.-6.*xw2/5.+0.5*xw3)
           ELSE
              phsmooth(i)=-one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-
     &                    xw3/30.)
              acsmooth(i)=-one/15.+8.*xw3/3.-3.*xw4+6.*xw3*xw2/5.-
     &                    xw4*xw2/6.
           ENDIF
           IF(xw.GE.two) THEN
              phsmooth(i)=one
              acsmooth(i)=one
           ENDIF
 30     CONTINUE

        IF(.NOT.inpmetal.AND..NOT.restart) THEN
           DO 40 p=1,nsph
              metals(p)=0.
 40        CONTINUE

           DO 50 p=nbodies-nstar+1,nbodies
              metals(p)=0.
 50        CONTINUE
        ENDIF

        IF(.NOT.inptform.AND..NOT.restart) THEN
           DO 60 p=nbodies-nstar+1,nbodies
              tform(p)=0.
 60        CONTINUE
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initbc
C
C
C***********************************************************************
C
C
C     Subroutine to set up the table of force and potential energy
C     correction terms due to the replicas arising from a particle 
C     located at the origin. Using symmetry, values are tabulated for
C     0 < x,y,z < 0.5 in the following arrays;
C
C          x-component of force correction -> bcfx
C          y-component of force correction -> bcfy
C          z-component of force correction -> bcfz
C          potential energy correction     -> bcphi
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER nx,ny,nz
        REAL deltar,x0,y0,z0,sx,sy,sz,phin
     
C=======================================================================

        deltar = 1./(2.*REAL(nmaxew))

        DO 30 nx = 0, nmaxew
           x0 = REAL(nx)*deltar
           DO 20 ny = 0, nx
              y0 = REAL(ny)*deltar
              DO 10 nz = 0, ny
                 z0 = REAL(nz)*deltar
                 CALL bcforce (x0,y0,z0,sx,sy,sz,phin)
C                     -------
                 bcfx(nx,ny,nz)  = sx
                 bcfy(nx,ny,nz)  = sy
                 bcfz(nx,ny,nz)  = sz
                 bcphi(nx,ny,nz) = phin
 10           CONTINUE
 20        CONTINUE
 30     CONTINUE
     
C-----------------------------------------------------------------------
C   Although its weight vanishes in subroutine ewaldmn, array(nmaxew+1)
C   is prepared for particles just on boundaries.
C-----------------------------------------------------------------------
     
        DO 50 ny = 0, nmaxew+1
           DO 40 nz = 0, ny
              bcfx(nmaxew+1,ny,nz)  = 0.0
              bcfy(nmaxew+1,ny,nz)  = 0.0
              bcfz(nmaxew+1,ny,nz)  = 0.0
              bcphi(nmaxew+1,ny,nz) = 0.0
 40        CONTINUE
 50     CONTINUE
     
C  Calculation of full arrays using symmetry.
C  ------------------------------------------
        DO 130 nx = 0, nmaxew+1
           DO 120 ny = 0, nmaxew+1
              DO 110 nz = 0, nmaxew+1
C  1. nx > nz > ny.
C  ----------------
                 IF(nx.GE.nz.AND.nz.GE.ny) THEN
                    bcfx(nx,ny,nz)  = bcfx(nx,nz,ny)
                    bcfy(nx,ny,nz)  = bcfz(nx,nz,ny)
                    bcfz(nx,ny,nz)  = bcfy(nx,nz,ny)
                    bcphi(nx,ny,nz) = bcphi(nx,nz,ny)
                    GO TO 110
                 ENDIF
C  2. ny > nx > nz.
C  ----------------
                 IF(ny.GE.nx.AND.nx.GE.nz) THEN
                    bcfx(nx,ny,nz)  = bcfy(ny,nx,nz)
                    bcfy(nx,ny,nz)  = bcfx(ny,nx,nz)
                    bcfz(nx,ny,nz)  = bcfz(ny,nx,nz)
                    bcphi(nx,ny,nz) = bcphi(ny,nx,nz)
                    GO TO 110
                 ENDIF
C  3. ny > nz > nx.
C  ----------------
                 IF(ny.GE.nz.AND.nz.GE.nx) THEN
                    bcfx(nx,ny,nz)  = bcfz(ny,nz,nx)
                    bcfy(nx,ny,nz)  = bcfx(ny,nz,nx)
                    bcfz(nx,ny,nz)  = bcfy(ny,nz,nx)
                    bcphi(nx,ny,nz) = bcphi(ny,nz,nx)
                    GO TO 110
                 ENDIF
C  4. nz > nx > ny.
C  ----------------
                 IF(nz.GE.nx.AND.nx.GE.ny) THEN
                    bcfx(nx,ny,nz)  = bcfy(nz,nx,ny)
                    bcfy(nx,ny,nz)  = bcfz(nz,nx,ny)
                    bcfz(nx,ny,nz)  = bcfx(nz,nx,ny)
                    bcphi(nx,ny,nz) = bcphi(nz,nx,ny)
                    GO TO 110
                 ENDIF
C  5. nz > ny > nx.
C  ----------------
                 IF(nz.GE.ny.AND.ny.GE.nx) THEN
                    bcfx(nx,ny,nz)  = bcfz(nz,ny,nx)
                    bcfy(nx,ny,nz)  = bcfy(nz,ny,nx)
                    bcfz(nx,ny,nz)  = bcfx(nz,ny,nx)
                    bcphi(nx,ny,nz) = bcphi(nz,ny,nx)
                 ENDIF
 110          CONTINUE
 120       CONTINUE
 130    CONTINUE
     
        RETURN
        END
C***********************************************************************
C
C
                SUBROUTINE bcforce(x,y,z,sx,sy,sz,phin)
C
C
C***********************************************************************
C
C     Correction force fields and potential energy due to the
C     periodic boundary condition are calculated according to
C     the Ewald summation method.
C        ( r_i - r_j(n) < (nreplica+0.6) L, h**2 < hmaxew2 )
C        (e.g., Sangster and Dixon, Adv.Phys. 25(1976)247.)
C
C     They are tabulated at (nmaxew+1)**3 grid points and in practice 
C     used by interpolation method in subroutine ewaldsum.
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER nx,ny,nz
        REAL x1,y1,z1,sx,sy,sz,phin,x0,x2,y0,y2,r2,z0,z2,ar2,ar,erbyr,
     &       ss0,hx,hx2,hy,hy2,h2,hz,hz2,ss,erfcc,x,y,z,r1,r3,reps,
     &       phi0,fx,fy,fz
     
C=======================================================================
     
        x1 = - x*2.0*pi
        y1 = - y*2.0*pi
        z1 = - z*2.0*pi
     
C-----------------------------------------------------------------------
C   A source is located at (x,y,z). The variables fx,fy,fz and phi are
C   force components and potential energy at (0,0,0) due to the source
C   at (x,y,z) + replicas. The variables sx,sy,sz and phin are force
C   components and potential energy at (0,0,0) due to the replica
C   sources only (without the direct term).
C-----------------------------------------------------------------------
     
        r2   = x*x + y*y + z*z

        IF(r2.LT.1.e-20) THEN
           sx = 0.
           sy = 0.
           sz = 0.
           phin = 0.
           RETURN
        ENDIF

        r1   = SQRT(r2)
        r3   = r1**3

        reps = r1/eps

        IF(reps.GE.2.0) THEN
           fx   = x/r3
           fy   = y/r3
           fz   = z/r3
           phi0 = 1./r1
        ELSE
          IF(reps.GE.1.0) THEN
             fx   = - 1.0/15.0+ reps**3*(8.0/3.0 + reps*(- 3.0 + 
     &              reps*(1.2 - reps/6.0)))
             fx   = fx/r3
             fy   = fx*y
             fz   = fx*z
             fx   = fx*x
             phi0 = 1.6 + reps*reps*(-4./3. + reps*(1. + reps*( - 0.3 + 
     &              reps/30.))) - 1./15./reps
             phi0 = phi0 / eps
          ELSE
             fx   = 4.0/3.0 + reps*reps*(-1.2 + reps*0.5)
             fx   = fx/(eps**3)
             fy   = fx*y
             fz   = fx*z
             fx   = fx*x
             phi0 = 1.4 + reps*reps*(-2./3. + reps*reps*(3.- reps)*0.1)
             phi0 = phi0 / eps
          ENDIF
        ENDIF

        sx   = - fx
        sy   = - fy
        sz   = - fz
        phin =  phi0
     
        DO 30  nx = -nreplica,nreplica
           x0 = x+REAL(nx)
           x2 = x0*x0
           DO 20 ny = -nreplica,nreplica
              y0 = y+REAL(ny)
              y2 = y0*y0
              r2 = x2 + y2
              IF( r2.gt.radew2 ) GO TO 20
              DO 10 nz = -nreplica,nreplica
                 z0 = z+REAL(nz)
                 z2 = z0*z0
                 r2 = x2 + y2 + z2
                 IF (r2.gt.radew2) GO TO 10
                 ar2 = alphaew*alphaew*r2
                 ar  = alphaew*SQRT(r2)
                 erbyr  = erfcc(ar)/SQRT(r2)
C                         -----
                 ss0 = erbyr/r2 + 2.0*alphaew*EXP(-ar2)/(sqrtpi*r2)
                 sx  = sx + x0*ss0
                 sy  = sy + y0*ss0
                 sz  = sz + z0*ss0
                 phin   = phin - erbyr
 10           CONTINUE
 20        CONTINUE
 30     CONTINUE
     
        DO 130  nx = -maxhew,maxhew
           hx  = REAL(nx)
           hx2 = hx*hx
           DO 120 ny = -maxhew,maxhew
              hy  = REAL(ny)
              hy2 = hy*hy
              h2  = hx2 + hy2
              IF( h2.gt.hmaxew2 ) GO TO 120
              DO 110 nz = -maxhew,maxhew
                 IF(nx.eq.0.AND.ny.eq.0.AND.nz.eq.0) GO TO 110
                 hz  = REAL(nz)
                 hz2 = hz*hz
                 h2  = hx2 + hy2 + hz2
                 IF (h2.gt.hmaxew2) GO TO 110
                 ss = 2.0*EXP(-h2*(pi/alphaew)**2)/h2*
     &                SIN(hx*x1+hy*y1+hz*z1)
                 sx = sx - hx*ss
                 sy = sy - hy*ss
                 sz = sz - hz*ss
                 phin = phin - EXP(-h2*(pi/alphaew)**2)/h2*
     &                  COS(hx*x1+hy*y1+hz*z1)/pi
 110          CONTINUE
 120       CONTINUE
 130    CONTINUE
     
        RETURN
        END
C***********************************************************************
C
C
                          FUNCTION erfcc(x)
C
C
C***********************************************************************
C
C
C      Function to compute the complementary error function erfc(x),
C      with fractional error everywhere less than 1.2x10^-7.  This
C      algorithm uses a Chebyshev fitting. (cf: Numerical Recipes
C      p.164).
C
C
C=======================================================================
     
C   Declaration of local variables.
C   -------------------------------
     
        REAL erfcc,x,z,t
     
C=======================================================================
     
        z = ABS(x)
        t = 1./(1.+0.5*z)
        erfcc = t*EXP( -z*z -1.26551223 + t*(1.00002368 + t*
     &               ( 0.37409196 + t*(0.09678418 + t*(-0.18628806 + t*
     &               (0.27886807 + t*(-1.13520398 + t*(1.48851587 + t*
     &               (-0.82215223 + t* 0.17087277   )))))))))
     
        IF (x.LT.0.0) erfcc = 2.0 - erfcc
     
        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initstep
C
C
C***********************************************************************
C
C
C     Subroutine to initialize parameters dealing with individual
C     particle time steps.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i

C=======================================================================

        DO 5 i=1,nbodies
           itimestp(i)=inittbin
           otimestp(i)=inittbin
 5      CONTINUE

        upbin=mintstep
        stime=0.0
        endstep=.TRUE.
        npactive=nbodies
        tsteppos=0.

        DO 55 i=1,npactive
           pactive(i)=i
 55     CONTINUE

        IF(usesph) THEN
           nsphact=nsph
        ELSE
           nsphact=0
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE zeroacc
C
C
C***********************************************************************
C
C
C     Subroutine to zero out acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i

C=======================================================================

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 10 i=1,npactive
           p=pactive(i)
           acc(p,1)=0.
           acc(p,2)=0.
           acc(p,3)=0.
 10     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE zeropot
C
C
C***********************************************************************
C
C
C     Subroutine to zero out potential.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i

C=======================================================================

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 10 i=1,npactive
           p=pactive(i)
           phi(p)=0.
           phiext(p)=0.
 10     CONTINUE

        IF(npactive.EQ.nbodies) esofttot=0.0

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initsph
C
C
C***********************************************************************
C
C
C     Subroutine to initialize properties of the SPH particles.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

        CALL maketsph
C            --------
        CALL inithsm
C            -------
        CALL neighbor('correct')
C            --------
        IF(sphinit) CALL initdens
C                        --------
        IF(sphinit.OR.ethinit.OR.readinas.OR.entinit) THEN
 
           IF(.NOT.uentropy) THEN

              CALL initeth
C                  -------
           ELSE

              CALL initent
C                  -------
           ENDIF
        ENDIF

        CALL initvet
C            -------
        CALL initdivv
C            --------
        IF(artfvisc.EQ.'bulk') THEN
           CALL accsphbv
C               --------
           IF(.NOT.isotherm) THEN

              IF(.NOT.uentropy) THEN

                 CALL ethdotbv('correct')
C                     --------
              ELSE

                 CALL entdotbv('correct')
C                     --------
              ENDIF
           ENDIF

        ELSE
           IF(artfvisc.EQ.'sphv') THEN
              CALL accsphcv
C                  --------
              IF(.NOT.isotherm) THEN

                 IF(.NOT.uentropy) THEN

                    CALL ethdotcv('correct')
C                        --------
                 ELSE

                    CALL entdotcv('correct')
C                        --------
                 ENDIF
              ENDIF
           ELSE
              CALL accsph
C                  ------
              IF(.NOT.isotherm) THEN

                 IF(.NOT.uentropy) THEN

                    CALL ethdot('correct')
C                        ------
                 ELSE

                    CALL entdot('correct')
C                        ------
                 ENDIF
              ENDIF

           ENDIF
        ENDIF

        IF(isotherm) CALL isomumax
C                         --------

        IF(radiate.AND.(.NOT.isotherm)) CALL ethrad
C                                            ------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE maketsph
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the tree structure for nearest 
C     neighbor searching.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

C-----------------------------------------------------------------------
C   Do not reinitialize tree if all particles are dissipative.
C-----------------------------------------------------------------------

        IF(nsph.NE.nbodies.OR.ttree.NE.tpos) THEN

           CALL setbox('sph ')
c..placeholder
C               ------

           CALL loadtree('sph ')
C               --------

           CALL cellnumb
C               --------

        ENDIF

C   Compute coordinates of edges of cells.
C   --------------------------------------
        CALL celledge
C            --------

C   Determine which cells contain ingroup particles.
C   ------------------------------------------------
        CALL groupcel
C            --------

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE setbox(option)
C
C
C***********************************************************************
C
C
C     Subroutine to adjust system box so that it contains all bodies.
C     The argument option indicates which particles to consider:
C     'all ' for all particles, 'sph ' for SPH particles, or 'coll'
C     for collisionless particles.  The local variable rebox indicates 
C     whether a resizing of the system box is to take place.  The 
C     variables posmin and posmax are the minimum and maximum 
C     coordinates of bodies in each dimension.  
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option
        INTEGER k,ismin,ismax,p
        LOGICAL rebox
        REAL posmin(ndim),posmax(ndim),posx(nbodsmax),posy(nbodsmax),
     &       posz(nbodsmax)

        EQUIVALENCE (posx(1),pos(1,1)),(posy(1),pos(1,2)),
     &              (posz(1),pos(1,3))

        SAVE rebox

        DATA rebox/.TRUE./

C=======================================================================

C   Determine minimum and maximum coordinates of bodies.
C   ----------------------------------------------------

        IF(boundary.EQ.'periodic') THEN

           rsize=pboxsize
           rmin(1)=-0.5*pboxsize
           rmin(2)=-0.5*pboxsize
           rmin(3)=-0.5*pboxsize

        ELSE

           IF(option.EQ.'all ') THEN
        
              posmin(1)=posx(ISMIN(nbodies,posx,1))
              posmin(2)=posy(ISMIN(nbodies,posy,1))
              posmin(3)=posz(ISMIN(nbodies,posz,1))
              posmax(1)=posx(ISMAX(nbodies,posx,1))
              posmax(2)=posy(ISMAX(nbodies,posy,1))
              posmax(3)=posz(ISMAX(nbodies,posz,1))

           ELSE

              IF(option.EQ.'sph ') THEN

                 posmin(1)=posx(ISMIN(nsph,posx,1))
                 posmin(2)=posy(ISMIN(nsph,posy,1))
                 posmin(3)=posz(ISMIN(nsph,posz,1))
                 posmax(1)=posx(ISMAX(nsph,posx,1))
                 posmax(2)=posy(ISMAX(nsph,posy,1))
                 posmax(3)=posz(ISMAX(nsph,posz,1))
 
              ELSE

                 DO 20 k=1,ndim

                    DO 10 p=nsph+1,nbodies
                       tempvect(p-nsph)=pos(p,k)
 10                 CONTINUE

                    posmin(k)=tempvect(ISMIN(nbodies-nsph,tempvect,1))
                    posmax(k)=tempvect(ISMAX(nbodies-nsph,tempvect,1))

 20              CONTINUE
   
              ENDIF
           ENDIF

C   Determine if a resizing is required.
C   ------------------------------------
           DO 50 k=1,ndim
              IF(rmin(k).GT.posmin(k).OR.rmin(k)+rsize.LT.posmax(k)) 
     &          rebox=.TRUE.
 50        CONTINUE

C   If a resizing is necessary, recompute rsize and rmin.
C   -----------------------------------------------------

           IF(rebox) THEN

              DO 70 k=1,ndim
                 IF(rsize.LT.posmax(k)-posmin(k)) rsize=posmax(k)-
     &              posmin(k)
 70           CONTINUE

              DO 80 k=1,ndim
                 rmin(k)=0.5*(posmin(k)+posmax(k))-0.5*rsize
 80           CONTINUE

           ENDIF

           rebox=.FALSE.

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE loadtree(option)
C
C
C***********************************************************************
C
C
C     Subroutine to insert the bodies into the tree.  The argument
C     option indicates if all bodies are lo be loaded ('all '), only
C     SPH particles are to be loaded ('sph '), or if only collisionless
C     particles are to be loaded ('coll').  The process is vectorized 
C     over active bodies.  Active bodies are those which are not yet in
C     place in the tree, as leaves.  The local variables pm1 and nindex
C     are used to convert back and forth between physical coordinates 
C     and subcell coordinates.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option
        INTEGER k,p,nindex(ndim),j,i,nbodlist,nclist,nclist2,
     &          nsubset,indcell,nsubbod1,nbodtemp,iupper,ilower,
     &          ncl2nsub
        REAL pm1(nsubcell,ndim)

        SAVE nindex,pm1

        DATA pm1/4*-1.,4*1.,2*-1.,2*1.,2*-1.,2*1.,-1.,1.,-1.,1.,
     &            -1.,1.,-1.,1./,nindex/4,2,1/

C=======================================================================

C   Deallocate old tree, compute coordinates of center of root cell.
C   ----------------------------------------------------------------
        incells=1
        root=nbodsmax+1

        DO 5 j=1,nsubcell
           subp(root,j)=0
 5      CONTINUE

        cellsize(root)=rsize

        DO 10 k=1,ndim
           pos(root,k)=rmin(k)+0.5*rsize
           bottom(root,k)=rmin(k)
 10     CONTINUE

C-----------------------------------------------------------------------
C   Place all bodies on active body list, having root as parent; place
C   root on active cell list.
C-----------------------------------------------------------------------

        IF(option.EQ.'all ') THEN
           ilower=1
           iupper=nbodies
        ELSE
           IF(option.EQ.'sph ') THEN
              ilower=1
              iupper=nsph
           ELSE
              ilower=nsph+1
              iupper=nbodies
           ENDIF
        ENDIF

        DO 20 i=ilower,iupper
           parent(i-ilower+1)=root
           bodlist(i-ilower+1)=i
 20     CONTINUE

        nbodlist=iupper-ilower+1
        celllist(1)=root
        nclist=1

C   Loop until no bodies are left active.
C   -------------------------------------

 200    CONTINUE

        IF(nclist.GT.0) THEN

C   Compute subindices for all active bodies.
C   -----------------------------------------
           DO 30 i=1,nbodlist
              subindex(i)=1
 30        CONTINUE

           DO 50 k=1,ndim
              DO 40 i=1,nbodlist
                 IF(pos(bodlist(i),k).GE.pos(parent(i),k)) 
     &                  subindex(i)=subindex(i)+nindex(k)
 40           CONTINUE
 50        CONTINUE

C   Compute number of bodies in each subcell.
C   -----------------------------------------
           DO 60 i=1,nbodlist
              subp(parent(i),subindex(i))=subp(parent(i),subindex(i))+1
 60        CONTINUE

C-----------------------------------------------------------------------
C   Open all subcells with more than one body, placing them on active 
C   cell list.
C-----------------------------------------------------------------------
           nclist2=0

           DO 110 j=1,nsubcell

              DO 70 i=1,nclist
                 asubp(i)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(nclist,asubp,1,1,isubset,nsubset)

              incells=incells+nsubset

              IF(incells.GT.ncells.OR.incells.GT.nbodsmax) 
     &           CALL terror(' overflow in loadtree')
C                     ------
              indcell=incells-nsubset+nbodsmax

              ncl2nsub=nclist2+nsubset

              IF(ncl2nsub.GT.nbodsmax.OR.ncl2nsub.GT.ncells)
     &           CALL terror(' nclist2 overflow in loadtree ')
C                     ------

              DO 90 k=1,nsubcell
                 DO 80 i=1,nsubset
                    subp(indcell+i,k)=0
 80              CONTINUE
 90           CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 100 i=1,nsubset
                 p=indcell+i
                 asubp(i)=celllist(isubset(i))
                 subp(asubp(i),j)=p
                 cellsize(p)=cellsize(asubp(i))*0.5
                 templist(nclist2+i)=p
                 pos(p,1)=pos(asubp(i),1)+pm1(j,1)*0.5*cellsize(p)
                 pos(p,2)=pos(asubp(i),2)+pm1(j,2)*0.5*cellsize(p)
                 pos(p,3)=pos(asubp(i),3)+pm1(j,3)*0.5*cellsize(p)
                 bottom(p,1)=pos(p,1)-0.5*cellsize(p)
                 bottom(p,2)=pos(p,2)-0.5*cellsize(p)
                 bottom(p,3)=pos(p,3)-0.5*cellsize(p)
 100          CONTINUE

              nclist2=nclist2+nsubset

 110       CONTINUE

           nclist=nclist2
        
           DO 120 i=1,nclist
              celllist(i)=templist(i)
 120       CONTINUE

C   Find all subcells with one body; add bodies to tree.
C   ----------------------------------------------------
           DO 130 i=1,nbodlist
              templist(i)=ncells*(subindex(i)-1)+(parent(i)-nbodsmax)
              asubp(i)=subpvect(templist(i))
 130       CONTINUE

           CALL WHENEQ(nbodlist,asubp,1,1,isubset,nsubbod1)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 140 i=1,nsubbod1
              subpvect(templist(isubset(i)))=bodlist(isubset(i))
 140       CONTINUE

C   Place bodies in cells with more than one body on active list.
C   ------------------------------------------------------------

           CALL WHENIGT(nbodlist,asubp,1,1,isubset,nbodtemp)

           nbodlist=nbodtemp

           DO 150 i=1,nbodlist
              parent(i)=asubp(isubset(i))
              templist(i)=bodlist(isubset(i))
 150       CONTINUE

           DO 160 i=1,nbodlist
              bodlist(i)=templist(i)
 160       CONTINUE

           GO TO 200
        
        ENDIF 

        DO 230 i=1,nbodies
           bottom(i,1)=pos(i,1)
           bottom(i,2)=pos(i,2)
           bottom(i,3)=pos(i,3)
 230    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE cellnumb
C
C
C***********************************************************************
C
C
C     Subroutine to compute the number of particles per cell for the
C     nearest neighbor search.  Vectorization is achieved by 
C     simultaneously processing all cells at the same level in the
C     hierarchy.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,nnodes,nsubc,nsubb

C=======================================================================

C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize cell properties.
C   ---------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           npercell(p)=0
 10     CONTINUE

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 20     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 30 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 40
              ENDIF
 30        CONTINUE

 40        CONTINUE

C   Compute properties of selected cells, looping over subcells.
C   ------------------------------------------------------------

           DO 100 j=1,nsubcell

              DO 50 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 50           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

              IF(nnodes.GT.ncells.OR.nnodes.GT.nbodsmax)
     &           CALL terror(' array overflow in cellnumb ')
C                     ------

              DO 60 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
 60           CONTINUE

              CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)               

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 70 i=1,nsubc
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+
     &              npercell(asubp(isubset(i)))
 70           CONTINUE

              CALL WHENILE(nnodes,asubp,1,nbodsmax,isubset,nsubb)

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nsubb
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+1
 80           CONTINUE

 100       CONTINUE

           fcell=lcell+1

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE celledge
C
C
C***********************************************************************
C
C
C     Subroutine to initialize coordinates of edges of cells for 
C     nearest neighbor searching.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER isgnsub(8,3),i,j,k,nclist,nclist2,nsubc,indlist

        SAVE isgnsub

        DATA isgnsub/4*0,4*1,0,0,1,1,0,0,1,1,0,1,0,1,0,1,0,1/

C=======================================================================

C   Initialize properties of root cell.
C   -----------------------------------

        DO 10 k=1,ndim
           bottom(root,k)=rmin(k)
 10     CONTINUE

C   Place root on list of cells.
C   ----------------------------
        celllist(1)=root
        nclist=1

C   Loop until all cells are processed.
C   -----------------------------------

 200    CONTINUE

        IF(nclist.GT.0) THEN

C   Select non-zero subcells from celllist.
C   ---------------------------------------
           nclist2=0

           DO 40 j=1,nsubcell

              DO 20 i=1,nclist
                 asubp(i)=subp(celllist(i),j)
 20           CONTINUE

              CALL WHENIGT(nclist,asubp,1,nbodsmax,isubset,nsubc)
  
              nclist2=nclist2+nsubc

              IF(nclist2.GT.ncells.OR.nclist2.GT.nbodsmax)
     &           CALL terror(' array overflow in celledge ')
C                     ------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nsubc
                 indlist=nclist2-nsubc+i
                 templist(indlist)=asubp(isubset(i))
                 parent(indlist)=celllist(isubset(i))
                 subindex(indlist)=j
 30           CONTINUE

 40        CONTINUE

           nclist=nclist2
        
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nclist
              celllist(i)=templist(i)
              bottom(celllist(i),1)=bottom(parent(i),1)+
     &                      isgnsub(subindex(i),1)*cellsize(celllist(i))    
              bottom(celllist(i),2)=bottom(parent(i),2)+
     &                      isgnsub(subindex(i),2)*cellsize(celllist(i))    
              bottom(celllist(i),3)=bottom(parent(i),3)+
     &                      isgnsub(subindex(i),3)*cellsize(celllist(i))    
 60        CONTINUE

           GO TO 200
        
        ENDIF 

        DO 70 i=1,nbodies
           bottom(i,1)=pos(i,1)
           bottom(i,2)=pos(i,2)
           bottom(i,3)=pos(i,3)
 70     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE groupcel
C
C
C***********************************************************************
C
C
C     Subroutine to locate cells containing ingroup particles, subject
C     to the constraint that the smoothing lengths of the particles
C     do not span too large a range.  Lists of such cells and the 
C     particles that lie within them are created. 
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER pstack,stack(nbodsmax),spointer,i,j,iblist,gpointer,
     &          gstack(nbodsmax),pgstack,numgroup,ismin,ismax
        LOGICAL testsubd
        REAL hmax,hmin

        EQUIVALENCE (stack(1),parent(1)),(gstack(1),asubp(1))

C=======================================================================

C   Locate cells containing ingroup particles.
C   ------------------------------------------
        stack(1)=root
        spointer=1
        ngroups=0
        iblist=0

 120    CONTINUE

        pstack=stack(spointer)
        spointer=spointer-1

C   If ingroup exceeded, subdivide cell; otherwise record cell.
C   -----------------------------------------------------------

        IF(pstack.LE.nbodsmax) THEN
           ngroups=ngroups+1
           groups(ngroups)=pstack
           iblist=iblist+1
           groupbod(iblist)=pstack
           pgroupb(ngroups)=iblist
        ELSE
           IF(npercell(pstack).GT.ingroup) THEN
              DO 130 j=1,nsubcell
                 IF(subp(pstack,j).GT.0) THEN
                    spointer=spointer+1
                    IF(spointer.GT.nbodsmax)
     &                 CALL terror(' array overflow in groupcel ')
C                           ------
                    stack(spointer)=subp(pstack,j)
                 ENDIF
 130          CONTINUE
           ELSE

C   Walk through tree, listing bodies within grouped cells.
C   -------------------------------------------------------

              gpointer=1
              gstack(1)=pstack
              numgroup=0

 150          CONTINUE

              pgstack=gstack(gpointer)
              gpointer=gpointer-1

              IF(pgstack.GT.nbodsmax) THEN

                 DO 160 j=1,nsubcell
                    IF(subp(pgstack,j).GT.0) THEN
                       gpointer=gpointer+1
                       IF(gpointer.GT.nbodsmax)
     &                    CALL terror(' array overflow in groupcel ')
C                              ------
                       gstack(gpointer)=subp(pgstack,j)
                    ENDIF
 160             CONTINUE

              ELSE

                 numgroup=numgroup+1
                 bodlist(numgroup)=pgstack

              ENDIF

              IF(gpointer.GT.0) GO TO 150

C   Subdivide if range in smoothing or softening lengths too large.
C   ---------------------------------------------------------------

              DO 170 i=1,numgroup
                 IF(bodlist(i).LE.nsph) THEN
                    tempvect(i)=hsmooth(bodlist(i))
                 ELSE
                    tempvect(i)=epsvect(bodlist(i))
                 ENDIF
 170          CONTINUE

              hmax=tempvect(ISMAX(numgroup,tempvect,1))
              hmin=tempvect(ISMIN(numgroup,tempvect,1))

              testsubd=cellsize(pstack).GT.5.*hmax.AND.hmax.GT.0.
              testsubd=testsubd.OR.hmax.GT.2.*hmin

              IF(testsubd) THEN

                 DO 180 j=1,nsubcell
                    IF(subp(pstack,j).GT.0) THEN
                       spointer=spointer+1
                       IF(spointer.GT.nbodsmax)
     &                    CALL terror(' array overflow in groupcel ')
C                              ------
                       stack(spointer)=subp(pstack,j)
                    ENDIF
 180             CONTINUE

              ELSE

                 ngroups=ngroups+1
                 groups(ngroups)=pstack
                 pgroupb(ngroups)=iblist+1

                 DO 190 i=1,numgroup
                    iblist=iblist+1
                    groupbod(iblist)=bodlist(i)
 190             CONTINUE

              ENDIF

           ENDIF
        ENDIF

        IF(spointer.GT.0) GO TO 120

        pgroupb(ngroups+1)=iblist+1

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE neighcol(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute number of near neighbors for all
C     collisionless particles.  The list of neighbors is returned 
C     from the subroutine findnear in the vector nearlist, which is 
C     equivalenced to the common array bodlist.  If argument pc is
C     'predict' then this subroutine performs a simple neighbor
C     search for a given set of softening lengths.  If pc is 'correct'
C     this subroutine will additionally adjust softening lengths so
C     that the number of neighbors is approximately nsvolume.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,pbody,npnear,np,inear,jsubset(nbodsmax),
     &          nearlist(nbodsmax),keepnear(nbodsmax),
     &          nearpb(nbodsmax),nnearbod
        LOGICAL savecrit
        REAL epsinv,testnn,dx,dy,dz

        EQUIVALENCE (nearlist(1),bodlist(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1)),(nearpb(1),parent(1))

C=======================================================================

C   Initialize neighbor diagnostics.
C   --------------------------------
        nstot=0
        nsmin=nbodies
        nsmax=0

C   Find nearest neighbors of grouped cells.
C   ----------------------------------------

        DO 100 p=1,ngroups

C   Find nearest neighbors.
C   -----------------------
           CALL findnear(p,npnear)
C               --------

           DO 70 np=pgroupb(p),pgroupb(p+1)-1

              pbody=groupbod(np)
              epsinv=1./epsvect(pbody)

              DO 20 i=1,npnear
                 dx=pos(pbody,1)-pos(nearlist(i),1)
                 dy=pos(pbody,2)-pos(nearlist(i),2)
                 dz=pos(pbody,3)-pos(nearlist(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 tempvect(i)=(dx**2+dy**2+dz**2)*epsinv*epsinv
                 IF(nearlist(i).EQ.pbody) tempvect(i)=5.0
 20           CONTINUE

C   Filter out neighbors further than two smoothing lengths.
C   --------------------------------------------------------
              CALL WHENFLT(npnear,tempvect,1,four,isubset,inear)

              testnn=ABS(REAL(inear-nsvolume)/REAL(nsvolume))
              savecrit=(.NOT.variabls).OR.(testnn.LE.nsvtol).OR.
     &                 (pc.EQ.'predict')

              IF(savecrit) THEN

                 nstot=nstot+inear
                 IF(inear.LT.nsmin) nsmin=inear
                 IF(inear.GT.nsmax) nsmax=inear
                 nnear(pbody)=inear

              ELSE

                 IF(inear.GT.nsvolume) THEN

                    DO 50 i=1,inear
                       nearpb(i)=nearlist(isubset(i))
 50                 CONTINUE

                    CALL reducee(pbody,inear,nnearbod)
C                        -------

                 ELSE

                    CALL enlargee(pbody,nnearbod)
C                        --------
                 ENDIF

              ENDIF

 70        CONTINUE

 100    CONTINUE

        IF(nbodies.NE.nsph) THEN
           nsavg=nstot/(nbodies-nsph)
        ELSE
           nsavg=0
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE inithsm
C
C
C***********************************************************************
C
C
C     Subroutine to initialize smoothing lengths and gravitational
C     softening lengths of SPH particles, allowing for a target number
C     of near neighbors.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j,nnearfix
        REAL avgdens,hsmthavg

C=======================================================================

        IF(sphinit) THEN

           IF(variablh.OR.consthsm.EQ.zero.OR.variablg.OR.epsgas.EQ.
     &        zero) THEN

              DO 10 i=1,nsph
                 avgdens=nsph/rsize**3
                 hsmooth(i)=0.01*(nsmooth/(4.*3.141592654*avgdens/
     &                      3.))**(one/3.)
 10           CONTINUE

              DO 40 j=1,20

                 CALL neighbor('predict')
C                     --------

                 DO 30 i=1,nsph
                    IF(nnear(i).EQ.0) THEN
                       nnearfix=1
                    ELSE
                       nnearfix=0
                    ENDIF
                    hsmooth(i)=0.5*hsmooth(i)*((REAL(nsmooth)/
     &                         REAL(nnear(i)+nnearfix))**(one/3.)+one)
 30              CONTINUE

 40           CONTINUE           

              IF(variablg) THEN
                 DO 45 i=1,nsph
                    epsvect(i)=hsmooth(i)
 45              CONTINUE
              ENDIF

              IF(((.NOT.variablh).AND.consthsm.EQ.zero).OR.
     &           ((.NOT.variablg).AND.epsgas.EQ.zero)) THEN

                 hsmthavg=0.

                 DO 50 i=1,nsph
                    hsmthavg=hsmthavg+hsmooth(i)
 50              CONTINUE

                 hsmthavg=hsmthavg/nsph

                 IF((.NOT.variablh).AND.consthsm.EQ.zero) THEN
                    DO 60 i=1,nsph
                       hsmooth(i)=hsmthavg
 60                 CONTINUE
                 ENDIF

                 IF((.NOT.variablg).AND.epsgas.EQ.zero) THEN
                    DO 65 i=1,nsph
                       epsvect(i)=hsmthavg
 65                 CONTINUE
                 ENDIF

              ENDIF

           ENDIF

           IF((.NOT.variablh).AND.consthsm.GT.zero) THEN

              DO 70 i=1,nsph
                 hsmooth(i)=consthsm
 70           CONTINUE

           ENDIF           

           IF((.NOT.variablg).AND.epsgas.GT.zero) THEN

              DO 80 i=1,nsph
                 epsvect(i)=epsgas
 80           CONTINUE

           ENDIF           

        ELSE

           IF(variablg) THEN

              DO 90 i=1,nsph
                 epsvect(i)=hsmooth(i)
 90           CONTINUE

           ELSE

              IF(epsgas.GT.zero) THEN

                 DO 100 i=1,nsph
                    epsvect(i)=epsgas
 100             CONTINUE

              ELSE

                 hsmthavg=0.

                 DO 110 i=1,nsph
                    hsmthavg=hsmthavg+hsmooth(i)
 110             CONTINUE

                 hsmthavg=hsmthavg/nsph

                 DO 120 i=1,nsph
                    epsvect(i)=hsmthavg
 120             CONTINUE

              ENDIF

           ENDIF

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initdens
C
C
C***********************************************************************
C
C
C     Subroutine to initialize smoothed estimate of the local density
C     for dissipative particles.  In addition, initdens smooths the
C     velocity field to provide relatively quiet initial conditions.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),iwsm,iwsm1
        REAL hsminv,wnorm,distnorm,dr2,dr2p,drw,wsm,wmass,dr2i,drw1,
     &       wsm1,wmass1,dx,dy,dz

        EQUIVALENCE (nb(1),bodlist(1))

C=======================================================================

C   Zero density, accumulators for velocity and thermal energy.
C   -----------------------------------------------------------
        DO 10 p=1,nsph
           rho(p)=0.
           veltpos(p,1)=vel(p,1)
           veltpos(p,2)=vel(p,2)
           veltpos(p,3)=vel(p,3)
           vel(p,1)=0.
           vel(p,2)=0.
           vel(p,3)=0.
 10     CONTINUE

C   Compute density of particles.
C   -----------------------------

        DO 30 p=1,nsph

           IF(symmetry.EQ.'hk') THEN
              hsminv=1.0/hsmooth(p)
              wnorm=piinv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

cc....        deldr2=4./ninterp=distance between interpolation points on kernel
cc....        deldr2i=1./deldr2
cc....        pnear(p)=pointer to the list of near bodies.
cc....        distnorm is 1/the physical distance between kernel interp pts.
cc....        wsmooth is the table of look-up values for the smoothing kernel

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 20 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 dr2=dx**2+dy**2+dz**2
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                 wmass=wnorm*wsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 wsm1=(1.-drw1)*wsmooth(iwsm1)+drw1*wsmooth(1+iwsm1)
                 wmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                  hsmooth(nb(i)))*wsm1
                 rho(p)=rho(p)+0.5*mass(nb(i))*(wmass+wmass1)
                 vel(p,1)=vel(p,1)+0.5*mass(nb(i))*(wmass+wmass1)*
     &                    veltpos(nb(i),1)
                 vel(p,2)=vel(p,2)+0.5*mass(nb(i))*(wmass+wmass1)*
     &                    veltpos(nb(i),2)
                 vel(p,3)=vel(p,3)+0.5*mass(nb(i))*(wmass+wmass1)*
     &                    veltpos(nb(i),3)
                 rho(nb(i))=rho(nb(i))+0.5*mass(p)*(wmass+wmass1)
                 vel(nb(i),1)=vel(nb(i),1)+0.5*mass(p)*(wmass+wmass1)*
     &                        veltpos(p,1)
                 vel(nb(i),2)=vel(nb(i),2)+0.5*mass(p)*(wmass+wmass1)*
     &                        veltpos(p,2)
                 vel(nb(i),3)=vel(nb(i),3)+0.5*mass(p)*(wmass+wmass1)*
     &                        veltpos(p,3)
 20           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 25 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx**2+dy**2+dz**2
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 wnorm=piinv*hsminv*hsminv*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                 wmass=wnorm*wsm
                 rho(p)=rho(p)+mass(nb(i))*wmass
                 vel(p,1)=vel(p,1)+mass(nb(i))*wmass*veltpos(nb(i),1)
                 vel(p,2)=vel(p,2)+mass(nb(i))*wmass*veltpos(nb(i),2)
                 vel(p,3)=vel(p,3)+mass(nb(i))*wmass*veltpos(nb(i),3)
                 rho(nb(i))=rho(nb(i))+mass(p)*wmass
                 vel(nb(i),1)=vel(nb(i),1)+mass(p)*wmass*veltpos(p,1)
                 vel(nb(i),2)=vel(nb(i),2)+mass(p)*wmass*veltpos(p,2)
                 vel(nb(i),3)=vel(nb(i),3)+mass(p)*wmass*veltpos(p,3)
 25           CONTINUE

           ENDIF

C   Include self-interaction term.
C   ------------------------------
           IF(symmetry.EQ.'be') THEN
              wnorm=piinv/hsmooth(p)**3
           ENDIF

           wmass=wnorm*mass(p)
           rho(p)=rho(p)+wmass
           vel(p,1)=vel(p,1)+wmass*veltpos(p,1)
           vel(p,2)=vel(p,2)+wmass*veltpos(p,2)
           vel(p,3)=vel(p,3)+wmass*veltpos(p,3)

 30     CONTINUE

        DO 40 p=1,nsph
           vel(p,1)=vel(p,1)/rho(p)
           vel(p,2)=vel(p,2)/rho(p)
           vel(p,3)=vel(p,3)/rho(p)
 40     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initeth
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the thermal energy.  Several options
C     are currently implemented.  If readinas is .TRUE., then the
C     thermal energy, u, is computed from the entropic function, a(s).
C     The two are related from the equation of state by: 
C
C               u = a(s) rho**(gamma-1) / (gamma-1).
C
C     If ethinit is .TRUE., then it is assumed that the thermal
C     energy has not been input and is to be computed using a Jeans
C     mass criterion.
C
C     If sphinit is .TRUE., then the thermal energy is smoothed
C     to provide relatively quiet initial conditions.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),iwsm,iwsm1
        REAL hsminv,wnorm,distnorm,dr2,dr2p,drw,wsm,wmass,dr2i,drw1,
     &       wsm1,wmass1,tcutoff,dx,dy,dz

        EQUIVALENCE (nb(1),bodlist(1))

C=======================================================================

        IF(ethinit) THEN

           DO 20 p=1,nsph
              tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(p)*nsmooth/8./ 
     &          (rho(p)/cosmof3))**(2./3.))-((epsvect(p)/  
     &          piinv)**2.)/3.)*(rho(p)/cosmof3) 
              IF(tcutoff.LT.mintemp) tcutoff=mintemp
              ethermal(p)=tcutoff/(meanmwt*mhboltz*gamma1)
 20        CONTINUE

        ENDIF

        IF(sphinit) THEN

C   Zero accumulator thermal energy.
C   --------------------------------
           DO 30 p=1,nsph
              ethold(p)=ethermal(p)
              ethermal(p)=0.
 30        CONTINUE

C   Compute density of particles.
C   -----------------------------

           DO 50 p=1,nsph

              IF(symmetry.EQ.'hk') THEN
                 hsminv=1.0/hsmooth(p)
                 wnorm=piinv*hsminv*hsminv*hsminv
                 distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 40 i=1,nnearlis(p)
                    nb(i)=nearbods(pnear(p)+i-1)
                    dx=pos(p,1)-pos(nb(i),1)
                    dy=pos(p,2)-pos(nb(i),2)
                    dz=pos(p,3)-pos(nb(i),3)
                    IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                    IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                    IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                    IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                    IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                    IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                    dr2=dx**2+dy**2+dz**2
                    dr2p=dr2*distnorm
                    iwsm=dr2p
                    IF(ninterp.LT.iwsm) iwsm=ninterp
                    drw=dr2p-iwsm
                    wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                    wmass=wnorm*wsm
                    dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                    iwsm1=dr2i
                    IF(ninterp.LT.iwsm1) iwsm1=ninterp
                    drw1=dr2i-iwsm1
                    wsm1=(1.-drw1)*wsmooth(iwsm1)+drw1*wsmooth(1+iwsm1)
                    wmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                     hsmooth(nb(i)))*wsm1
                    ethermal(p)=ethermal(p)+0.5*mass(nb(i))*(wmass+
     &                          wmass1)*ethold(nb(i))
                    ethermal(nb(i))=ethermal(nb(i))+0.5*mass(p)*
     &                              (wmass+wmass1)*ethold(p)
 40              CONTINUE

              ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 45 i=1,nnearlis(p)
                    nb(i)=nearbods(pnear(p)+i-1)
                    dx=pos(p,1)-pos(nb(i),1)
                    dy=pos(p,2)-pos(nb(i),2)
                    dz=pos(p,3)-pos(nb(i),3)
                    IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                    IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                    IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                    IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                    IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                    IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                    dr2=dx**2+dy**2+dz**2
                    hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                    wnorm=piinv*hsminv*hsminv*hsminv
                    distnorm=hsminv**2*deldr2i
                    dr2p=dr2*distnorm   
                    iwsm=dr2p
                    IF(ninterp.LT.iwsm) iwsm=ninterp
                    drw=dr2p-iwsm
                    wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                    wmass=wnorm*wsm
                    ethermal(p)=ethermal(p)+mass(nb(i))*wmass*
     &                          ethold(nb(i))
                    ethermal(nb(i))=ethermal(nb(i))+mass(p)*
     &                              wmass*ethold(p)
 45              CONTINUE

              ENDIF

C   Include self-interaction term.
C   ------------------------------
              IF(symmetry.EQ.'be') THEN
                 wnorm=piinv/hsmooth(p)**3
              ENDIF

              wmass=wnorm*mass(p)
              ethermal(p)=ethermal(p)+wmass*ethold(p)

 50        CONTINUE

           DO 60 p=1,nsph
              ethermal(p)=ethermal(p)/rho(p)
 60        CONTINUE

        ENDIF

        IF(readinas.AND.(.NOT.ethinit)) THEN

           IF(.NOT.isotherm) THEN

              DO 70 p=1,nsph
                 ethermal(p)=ethermal(p)*rho(p)**gamma1/gamma1
 70           CONTINUE

           ELSE

              CALL terror(' input error in initeth -- gamma = 1 ')
C                  ------

           ENDIF

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initvet
C
C
C***********************************************************************
C
C
C     Subroutine to initialize thermodynamic variables and velocity 
C     which is time-syncrhonized with the position.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k
        REAL rhot

C=======================================================================

C### LH 050990

        IF(cosmo.AND.comove) THEN
           DO 20 k=1,ndim
              DO 10 p=1,nsph
                 veltpos(p,k)=cosmofac*vel(p,k)
 10           CONTINUE
 20        CONTINUE
        ELSE
           DO 40 k=1,ndim
              DO 30 p=1,nsph
                 veltpos(p,k)=vel(p,k)
 30           CONTINUE
 40        CONTINUE
        ENDIF

C### LH 050990

        IF(.NOT.isotherm) THEN

           IF(.NOT.uentropy) THEN

              DO 50 p=1,nsph
                 csound(p)=SQRT(gamma*gamma1*ethermal(p))
                 dethdt(p)=0.
                 derad(p)=0.
 50           CONTINUE

           ELSE

              DO 60 p=1,nsph
                 rhot=rho(p)*cosmof3
                 csound(p)=SQRT(gamma*rhot**gamma1*entropy(p))
                 dentdt(p)=0.
                 derad(p)=0.
 60           CONTINUE

           ENDIF

        ELSE
           DO 70 p=1,nsph
              csound(p)=SQRT(ethermal(p))
 70        CONTINUE
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initdivv
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the divergence and magnitude of the 
C     curl of the velocity field for SPH particles.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),iwsm,iwsm1
        REAL hsminv,distnorm,dr2,dr2p,drw,dr2i,drw1,dwnorm,vdotdr,dwsm,
     &       dwmass,dwsm1,dwmass1,dx,dy,dz,dvx,dvy,dvz,curlvx(nsphmax),
     &       curlvy(nbodsmax),curlvz(nworkvec),dwmp,dwmnbi

        EQUIVALENCE (nb(1),bodlist(1)),(curlvx(1),hsmcurlv(1)),
     &              (curlvy(1),tempvect(1)),(curlvz(1),workvect(1))

C=======================================================================

C   Zero the velocity divergence.
C   -----------------------------
        DO 10 p=1,nsph
           hsmdivv(p)=0.
           curlvx(p)=0.
           curlvy(p)=0.
           curlvz(p)=0.
 10     CONTINUE

C   Compute div v of particles.
C   ---------------------------

        DO 30 p=1,nsph

           IF(symmetry.EQ.'hk') THEN
              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 20 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 dvx=veltpos(p,1)-veltpos(nb(i),1)+dx*cosmohub
                 dvy=veltpos(p,2)-veltpos(nb(i),2)+dy*cosmohub
                 dvz=veltpos(p,3)-veltpos(nb(i),3)+dz*cosmohub
                 vdotdr=dvx*dx+dvy*dy+dvz*dz
                 dwmp=mass(p)*0.5*(dwmass+dwmass1)
                 dwmnbi=mass(nb(i))*0.5*(dwmass+dwmass1)
                 hsmdivv(p)=hsmdivv(p)-dwmnbi*vdotdr
                 curlvx(p)=curlvx(p)+dwmnbi*(dz*dvy-dy*dvz)
                 curlvy(p)=curlvy(p)+dwmnbi*(dx*dvz-dz*dvx)
                 curlvz(p)=curlvz(p)+dwmnbi*(dy*dvx-dx*dvy)
                 hsmdivv(nb(i))=hsmdivv(nb(i))-dwmp*vdotdr
                 curlvx(nb(i))=curlvx(nb(i))+dwmp*(dz*dvy-dy*dvz)
                 curlvy(nb(i))=curlvy(nb(i))+dwmp*(dx*dvz-dz*dvx)
                 curlvz(nb(i))=curlvz(nb(i))+dwmp*(dy*dvx-dx*dvy)
 20           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 25 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dvx=veltpos(p,1)-veltpos(nb(i),1)+dx*cosmohub
                 dvy=veltpos(p,2)-veltpos(nb(i),2)+dy*cosmohub
                 dvz=veltpos(p,3)-veltpos(nb(i),3)+dz*cosmohub
                 vdotdr=dvx*dx+dvy*dy+dvz*dz
                 dwmp=mass(p)*dwmass
                 dwmnbi=mass(nb(i))*dwmass
                 hsmdivv(p)=hsmdivv(p)-dwmnbi*vdotdr
                 curlvx(p)=curlvx(p)+dwmnbi*(dz*dvy-dy*dvz)
                 curlvy(p)=curlvy(p)+dwmnbi*(dx*dvz-dz*dvx)
                 curlvz(p)=curlvz(p)+dwmnbi*(dy*dvx-dx*dvy)
                 hsmdivv(nb(i))=hsmdivv(nb(i))-dwmp*vdotdr
                 curlvx(nb(i))=curlvx(nb(i))+dwmp*(dz*dvy-dy*dvz)
                 curlvy(nb(i))=curlvy(nb(i))+dwmp*(dx*dvz-dz*dvx)
                 curlvz(nb(i))=curlvz(nb(i))+dwmp*(dy*dvx-dx*dvy)
 25           CONTINUE

           ENDIF

 30     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 p=1,nsph
           hsmdivv(p)=hsmooth(p)*hsmdivv(p)/rho(p)
           hsmcurlv(p)=hsmooth(p)*SQRT(curlvx(p)**2+curlvy(p)**2+
     &                 curlvz(p)**2)/rho(p)
 40     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE isomumax
C
C
C***********************************************************************
C
C
C     Subroutine to compute mumaxdvh if an isothermal equation of
C     state is specified.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),imumax,ismax
        REAL dx,dy,dz,dr2,muij(nbodsmax),vdotdr,tmuij,hsmavg,formp,
     &       formnbi,abshdivv,hsmdvvni,hsmdvvp

        EQUIVALENCE (nb(1),bodlist(1)),(muij(1),tempvect(1))

C=======================================================================

        DO 40 p=1,nsph

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

           IF(artfvisc.EQ.'sph '.OR.artfvisc.EQ.' sph') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 10 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 muij(i)=ABS(muij(i))
 10           CONTINUE

           ENDIF

           IF(artfvisc.EQ.'sphv') THEN

              abshdivv=ABS(hsmdivv(p))
              formp=abshdivv/(abshdivv+hsmcurlv(p)+0.0001*csound(p))

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 20 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/(dr2+epssph*
     &                   hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 muij(i)=ABS(muij(i))
 20           CONTINUE

           ENDIF

           IF(artfvisc.EQ.'bulk') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    muij(i)=hsmdvvni
                 ELSE
                    muij(i)=hsmdvvp
                 ENDIF
 30           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0) THEN
              imumax=ISMAX(nnearlis(p),muij,1)
              tmuij=muij(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 tmuij=muij(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 35           CONTINUE

           ENDIF

 40     CONTINUE

        teth=tnow

        RETURN
        END
C***********************************************************************
C
C
                            SUBROUTINE ethrad
C
C
C***********************************************************************
C
C
C     Subroutine to compute source and sink terms in the energy
C     equation resulting from radiative heating and cooling of the
C     gas.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p
        REAL temp,cool,templog,deradold,cutoff,tempdiff,tenf3,
     &       tenf3exp,compcool,dt,aconst,dtaconst,tcutoff,tcutoffc,
     &       compdiff,cutoffc

C=======================================================================

        dt=teth-trad
        dtaconst=dt
        IF(dt.EQ.0.) dtaconst=dtime/(2.*inittbin)

        DO 10 p=1,nsph
           deradold=derad(p)
           temp=meanmwt*mhboltz*csound(p)**2/gamma
           tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(p)*nsmooth/8./
     &          (rho(p)/cosmof3))**(2./3.))-((epsvect(p)/
     &          piinv)**2.)/3.)*(rho(p)/cosmof3)
           IF(tcutoff.LT.mintemp) tcutoff=mintemp
           templog=LOG10(temp)
           tempdiff= 30.*(1.-temp/tcutoff)
           IF(tempdiff.LT.0.)tempdiff=0.
           cutoff=exp(-tempdiff *tempdiff)

           tcutoffc=tcutoff
           IF(tcutoff.LT.10000.) tcutoffc=10000.
           compdiff= 30.*(1.-temp/tcutoffc)
           IF(compdiff.LE.0.)compdiff=0.
           cutoffc=exp(-compdiff *compdiff)

           tenf3exp=6.2-templog
           IF(templog-6.2.GT.zero) tenf3exp=zero
           tenf3=10.**(-2.*tenf3exp**4-1.7)
           cool=10.**(-4.43-0.273*(4.-templog)**2)+
     &          10.**(-0.10-1.880*(5.23-templog)**4)+tenf3
           compcool=comptmh*temp*((1.+redshift)**4)
           derad(p)=aheatmh*fhydrogn+bheatmh2*rho(p)*fhydrog2/
     &              cosmof3-ccoolmh2*rho(p)*cutoff*cool*
     &              fhydrog2/cosmof3-compcool*cutoffc

           IF(.NOT.uentropy) THEN
              aconst=slowcool*csound(p)**2/(gamma*gamma1*dtaconst)+
     &               dethdt(p)
           ELSE
              aconst=slowcool*csound(p)**2/(gamma*gamma1*dtaconst)+
     &               rho(p)**gamma1*dentdt(p)/gamma1
           ENDIF

           IF(aconst.LE.zero) aconst=tiny
           aconst=derad(p)/aconst
           derad(p)=derad(p)/SQRT(1.+aconst*aconst)
           eradiate=eradiate+0.5*dt*mass(p)*(derad(p)+deradold)
 10     CONTINUE

        trad=teth

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initeps
C
C
C***********************************************************************
C
C
C     Subroutine to initialize gravitational softening lengths for
C     collisionless particles, allowing for a target number of near 
C     neighbors.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,j,nnearfix
        REAL avgdens,epsavg

C=======================================================================

        IF(variabls.OR.eps.EQ.zero) THEN

           DO 10 i=nsph+1,nbodies
              avgdens=(nbodies-nsph)/rsize**3
              epsvect(i)=0.01*(nsvolume/(4.*3.141592654*avgdens/3.))**
     &                   (one/3.)
 10        CONTINUE

           DO 40 j=1,20

              CALL neighcol('predict')
C                  --------

              DO 30 i=nsph+1,nbodies
                 IF(nnear(i).EQ.0) THEN
                    nnearfix=1
                 ELSE
                    nnearfix=0
                 ENDIF
                 epsvect(i)=0.5*epsvect(i)*((REAL(nsvolume)/
     &                      REAL(nnear(i)+nnearfix))**(one/3.)+one)
 30           CONTINUE

 40        CONTINUE           

           IF(.NOT.variabls) THEN

              epsavg=0.

              DO 50 i=nsph+1,nbodies
                 epsavg=epsavg+epsvect(i)
 50           CONTINUE

              IF(nsph.NE.nbodies) epsavg=epsavg/(nbodies-nsph)

              DO 60 i=nsph+1,nbodies
                 epsvect(i)=epsavg
 60           CONTINUE

           ENDIF

        ELSE

           DO 70 i=nsph+1,nbodies
              epsvect(i)=eps
 70        CONTINUE

        ENDIF           

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initpos
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the positions of the bodies for the
C     initial timestep, dtime/initbin.  The local variable p is a
C     pointer to the bodies.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k
        REAL dt2,acceff,dt

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        dt=dtime/inittbin
        dt2=dt**2

        IF(boundary.EQ.'vacuum') THEN

           IF((.NOT.friction).OR.friclucy) THEN

              DO 200 k=1,ndim
                 DO 100 p=1,nbodies
                    pos(p,k)=pos(p,k)+acc(p,k)*dt2/8.
 100             CONTINUE
 200          CONTINUE

           ELSE

              DO 250 k=1,ndim
                 DO 230 p=1,nbodies
                    IF(p.LE.nsph) THEN
                       acceff=cfrict/dt
                    ELSE
                       acceff=zero
                    ENDIF
                    pos(p,k)=pos(p,k)+(acc(p,k)-(acceff*vel(p,k)))*
     &                       dt2/8.
 230             CONTINUE
 250          CONTINUE

           ENDIF

        ELSE

           IF(.NOT.comove) THEN

              DO 400 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 300 p=1,nbodies
                    pos(p,k)=pos(p,k)+acc(p,k)*dt2/8.
                    IF(pos(p,k).GE.hboxsize) pos(p,k)=pos(p,k)-pboxsize
                    IF(pos(p,k).LT.-hboxsize) pos(p,k)=pos(p,k)+pboxsize
 300             CONTINUE
 400          CONTINUE

           ELSE

              DO 600 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 500 p=1,nbodies
                    pos(p,k)=pos(p,k)+acc(p,k)*dt2/(8.*cosmof3)-
     &                       vel(p,k)*dt2*hubble/4.
                    IF(pos(p,k).GE.hboxsize) pos(p,k)=pos(p,k)-pboxsize
                    IF(pos(p,k).LT.-hboxsize) pos(p,k)=pos(p,k)+pboxsize
 500             CONTINUE
 600          CONTINUE

           ENDIF

        ENDIF

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE stepsys(n)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the state of the system by one large
C     timestep.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n

C=======================================================================


 50     CONTINUE

           CALL timestep
C               --------

C   Update positions.
C   -----------------
           CALL steppos
C               -------

           IF(npactive.GT.0) THEN

C   Obtain extrapolated estimate of velocities.
C   -------------------------------------------
              IF(usesph.AND.nsphact.GT.0) CALL vextrap
C                                              -------
C   Zero out acceleration, compute gravitational acceleration.
C   ----------------------------------------------------------
              IF(.NOT.endstep) THEN

                 CALL zeroacc
C                     -------
                 CALL gravity('acc ')
C                     -------
              ENDIF

C   Advance thermal energy, compute SPH forces.
C   -------------------------------------------
              IF(usesph.AND.nsphact.GT.0) CALL stepsph
C                                              -------
           ENDIF
 
           IF(.NOT.endstep) THEN

              CALL stepvel
C                  -------
              IF(sphfiltr.AND.usesph.AND.nsphact.GT.0) CALL vfilter
C                                                           -------
              GO TO 50

           ENDIF

C   Output system state.
C   --------------------
        CALL outstate(n)
C            --------

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE timestep
C
C
C***********************************************************************
C
C
C     Subroutine to find appropriate time step for each particle
C     for the variable time step scheme.  Also finds the smallest
C     timestep and determines which particles are to be moved during
C     the next time step
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,p,ttemp,ideal,nrem,npacts,npactl,tmin,upbint,ismax,
     &          nchange,isrchigt,idele,nsubset,nsubset1
        LOGICAL testcrit
        REAL epart,vt,at,c1,tflag,epartp,abshdivv

C=======================================================================

        IF(ektot.LT.-eptot/6.) THEN
           epart= -eptot/(6.*mtot)
        ELSE
           epart=ektot/mtot
        ENDIF

        CALL WHENIGT(nbodies,itimestp,1,upbin,pactive,npactive)

        nsphact=ISRCHIGT(npactive,pactive,1,nsph)-1
        IF(nsphact.LT.0) nsphact=0

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 6 i=1,nsphact
           p=pactive(i)
           abshdivv=hsmdivv(p)
           abshdivv=ABS(abshdivv)
           templist(i)=dtime*(csound(p)+1.2*(alpha*csound(p)+beta*
     &                 mumaxdvh(p))+abshdivv)/(hsmooth(p)*cosmofac*
     &                 courant)+1.
           mumaxdvh(p)=0.
 6      CONTINUE

        DO 7 i=nsphact+1,npactive
           templist(i)=0
 7      CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 10 i=1,npactive
           p=pactive(i)
           vt=vel(p,1)*vel(p,1)+vel(p,2)*vel(p,2)+vel(p,3)*vel(p,3)
           at=acc(p,1)*acc(p,1)+acc(p,2)*acc(p,2)+acc(p,3)*acc(p,3)
           epartp=0.5*vt+phi(p)+phiext(p)
           epartp=ABS(epartp)
           IF(cosmo.AND.(.NOT.comove)) THEN
              epartp=MAX(epart,epartp)
           ENDIF
           IF(epartp.NE.0.0) THEN
              idele=dtime*SQRT(at*vt)/(etol*epartp)+1.
           ELSE
              idele=1
           ENDIF
           IF(templist(i).LT.idele) templist(i)=idele
           testcrit=itimestp(p).LT.templist(i).AND.itimestp(p)
     &              .LT.mintstep
           IF(testcrit) THEN
              tflag=two
           ELSE
              tflag=zero
           ENDIF
           testcrit=(itimestp(p).NE.1).AND.(itimestp(p).GT.2*upbin)
     &              .AND.(itimestp(p)/2.GE.templist(i))
           IF(testcrit) THEN
              tempvect(i)=four
           ELSE
              tempvect(i)=zero
           ENDIF
           tempvect(i)=tempvect(i)+tflag
 10     CONTINUE

        CALL WHENFGT(npactive,tempvect,1,three,isubset,nchange)

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 15 i=1,nchange
           p=pactive(isubset(i))
           itimestp(p)=itimestp(p)/2
           tempvect(isubset(i))=0.
 15     CONTINUE

        CALL WHENFGT(npactive,tempvect,1,one,isubset,nchange)

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 16 i=1,nchange
           p=pactive(isubset(i))
           c1=templist(isubset(i))
           ideal=LOG(c1)/log2 + 1.
           ideal=2**ideal + .5
           IF(mintstep.LT.ideal) THEN
              itimestp(p)=mintstep
           ELSE
              itimestp(p)=ideal
           ENDIF
 16     CONTINUE

        DO 17 p=1,nbodies
           tempvect(p)=itimestp(p)
 17     CONTINUE

        tmin=itimestp(ISMAX(nbodies,tempvect,1))

        stime=stime + 1./(tmin*2.)
        upbin=tmin
        ttemp=2.*tmin*stime+.5

        IF(ttemp.GE.2*tmin) THEN

           tsteppos=dtime*(1.-stime+1./(tmin*2.))
           upbin=0
           stime=0.0
           npactive=nbodies
           endstep=.TRUE.

           DO 55 i=1,npactive
              pactive(i)=i
 55        CONTINUE

           nsphact=nsph

           RETURN
   
        ENDIF

 20     CONTINUE

        IF(upbin.GT.1.AND.MOD(ttemp,2).EQ.0) THEN
           upbin=upbin/2
           ttemp=ttemp/2
           GO TO 20
        ENDIF

        CALL WHENEQ(nbodies,itimestp,1,upbin,pactive,npactive)

        IF(npactive.GT.0) THEN

        upbint=upbin

 23        CONTINUE

           nrem=MOD(npactive,ntvector)
           nrem=ntvector-nrem
           npacts=MOD(nrem,ntvector)

           IF(npacts.NE.0.AND.upbint.NE.1) THEN
              npactl=npactive
              upbint=upbint/2
              IF(upbint.EQ.0) CALL terror('error in timestep')
C                                  ------

              CALL WHENEQ(nbodies,itimestp,1,upbint,bodlist,nchange)

              IF(nchange.LT.npacts) npacts=nchange
              npactive=npactive+npacts

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,npacts
                 p=bodlist(i)
                 pactive(i+npactl)=p
                 itimestp(p)=upbin
 30           CONTINUE
              GO TO 23
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,npactive
              p=pactive(i)
              isubset(i)=itimestp(p)/otimestp(p)
              otimestp(p)=itimestp(p)
 40        CONTINUE

           CALL WHENEQ(npactive,isubset,1,0,templist,nchange)

           IF((.NOT.comove).OR.(.NOT.cosmo)) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 45 i=1,nchange
                 p=pactive(templist(i))
                 c1=3./32.
                 c1=c1*dtime*dtime/(itimestp(p)*itimestp(p))
                 pos(p,1)=pos(p,1)+c1*acc(p,1)
                 pos(p,2)=pos(p,2)+c1*acc(p,2)
                 pos(p,3)=pos(p,3)+c1*acc(p,3)
                 IF(pos(p,1).GE.hboxsize.AND.bwrap) pos(p,1)=pos(p,1)-
     &                                                  pboxsize
                 IF(pos(p,1).LT.-hboxsize.AND.bwrap) pos(p,1)=pos(p,1)+
     &                                                  pboxsize
                 IF(pos(p,2).GE.hboxsize.AND.bwrap) pos(p,2)=pos(p,2)-
     &                                                  pboxsize
                 IF(pos(p,2).LT.-hboxsize.AND.bwrap) pos(p,2)=pos(p,2)+
     &                                                  pboxsize
                 IF(pos(p,3).GE.hboxsize.AND.bwrap) pos(p,3)=pos(p,3)-
     &                                                  pboxsize
                 IF(pos(p,3).LT.-hboxsize.AND.bwrap) pos(p,3)=pos(p,3)+
     &                                                  pboxsize
 45           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 47 i=1,nchange
                 p=pactive(templist(i))
                 c1=3./32.
                 c1=c1*dtime*dtime/(itimestp(p)*itimestp(p))
                 pos(p,1)=pos(p,1)+c1*acc(p,1)/cosmof3-2.*c1*
     &                    hubble*vel(p,1)
                 pos(p,2)=pos(p,2)+c1*acc(p,2)/cosmof3-2.*c1*
     &                    hubble*vel(p,2)
                 pos(p,3)=pos(p,3)+c1*acc(p,3)/cosmof3-2.*c1*
     &                    hubble*vel(p,3)
                 IF(pos(p,1).GE.hboxsize.AND.bwrap) pos(p,1)=pos(p,1)-
     &                                                  pboxsize
                 IF(pos(p,1).LT.-hboxsize.AND.bwrap) pos(p,1)=pos(p,1)+
     &                                                  pboxsize
                 IF(pos(p,2).GE.hboxsize.AND.bwrap) pos(p,2)=pos(p,2)-
     &                                                  pboxsize
                 IF(pos(p,2).LT.-hboxsize.AND.bwrap) pos(p,2)=pos(p,2)+
     &                                                  pboxsize
                 IF(pos(p,3).GE.hboxsize.AND.bwrap) pos(p,3)=pos(p,3)-
     &                                                  pboxsize
                 IF(pos(p,3).LT.-hboxsize.AND.bwrap) pos(p,3)=pos(p,3)+
     &                                                  pboxsize
 47           CONTINUE

           ENDIF

           CALL WHENIGT(npactive,isubset,1,1,templist,nchange)

           IF((.NOT.comove).OR.(.NOT.cosmo)) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 50 i=1,nchange
                 p=pactive(templist(i))
                 c1=-1./8.*(1.-1./isubset(templist(i)))*
     &            (1.+1./isubset(templist(i)))*(isubset(templist(i))**2)
                 c1=c1*dtime*dtime/(itimestp(p)*itimestp(p))
                 pos(p,1)=pos(p,1)+c1*acc(p,1)
                 pos(p,2)=pos(p,2)+c1*acc(p,2)
                 pos(p,3)=pos(p,3)+c1*acc(p,3)
                 IF(pos(p,1).GE.hboxsize.AND.bwrap) pos(p,1)=pos(p,1)-
     &                                                  pboxsize
                 IF(pos(p,1).LT.-hboxsize.AND.bwrap) pos(p,1)=pos(p,1)+
     &                                                  pboxsize
                 IF(pos(p,2).GE.hboxsize.AND.bwrap) pos(p,2)=pos(p,2)-
     &                                                  pboxsize
                 IF(pos(p,2).LT.-hboxsize.AND.bwrap) pos(p,2)=pos(p,2)+
     &                                                  pboxsize
                 IF(pos(p,3).GE.hboxsize.AND.bwrap) pos(p,3)=pos(p,3)-
     &                                                  pboxsize
                 IF(pos(p,3).LT.-hboxsize.AND.bwrap) pos(p,3)=pos(p,3)+
     &                                                  pboxsize
 50           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 52 i=1,nchange
                 p=pactive(templist(i))
                 c1=-1./8.*(1.-1./isubset(templist(i)))*
     &            (1.+1./isubset(templist(i)))*(isubset(templist(i))**2)
                 c1=c1*dtime*dtime/(itimestp(p)*itimestp(p))
                 p=pactive(templist(i))
                 pos(p,1)=pos(p,1)+c1*acc(p,1)/cosmof3-2.*c1*
     &                    hubble*vel(p,1)
                 pos(p,2)=pos(p,2)+c1*acc(p,2)/cosmof3-2.*c1*
     &                    hubble*vel(p,2)
                 pos(p,3)=pos(p,3)+c1*acc(p,3)/cosmof3-2.*c1*
     &                    hubble*vel(p,3)
                 IF(pos(p,1).GE.hboxsize.AND.bwrap) pos(p,1)=pos(p,1)-
     &                                                  pboxsize
                 IF(pos(p,1).LT.-hboxsize.AND.bwrap) pos(p,1)=pos(p,1)+
     &                                                  pboxsize
                 IF(pos(p,2).GE.hboxsize.AND.bwrap) pos(p,2)=pos(p,2)-
     &                                                  pboxsize
                 IF(pos(p,2).LT.-hboxsize.AND.bwrap) pos(p,2)=pos(p,2)+
     &                                                  pboxsize
                 IF(pos(p,3).GE.hboxsize.AND.bwrap) pos(p,3)=pos(p,3)-
     &                                                  pboxsize
                 IF(pos(p,3).LT.-hboxsize.AND.bwrap) pos(p,3)=pos(p,3)+
     &                                                  pboxsize
 52           CONTINUE

           ENDIF

        ENDIF

        tsteppos=dtime/(2.*tmin)
        endstep=.FALSE.

        DO 60 i=1,npactive
           templist(i)=pactive(i)
 60     CONTINUE

        CALL WHENILE(npactive,templist,1,nsph,isubset,nsubset)

        DO 70 i=1,nsubset
           pactive(i)=templist(isubset(i))
 70     CONTINUE

        CALL WHENIGT(npactive,templist,1,nsph,isubset,nsubset1)

        DO 80 i=1,nsubset1
           pactive(nsubset+i)=templist(isubset(i))
 80     CONTINUE

        nsphact=nsubset

        RETURN

        END
C***********************************************************************
C
C
                          SUBROUTINE steppos
C
C
C***********************************************************************
C
C
C     Subroutine to advance the positions of the bodies for a
C     timestep dt.  The local variable p is a pointer to the bodies.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k
        REAL acceff

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------

        IF((.NOT.friction).OR.(.NOT.friclucy)) THEN

              DO 200 k=1,ndim
                 DO 100 p=1,nbodies
                    pos(p,k)=pos(p,k)+vel(p,k)*tsteppos
 100             CONTINUE
 200          CONTINUE

        ELSE

           DO 600 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 500 p=1,nbodies
                 IF(p.LE.nsph) THEN
                    acceff=cfrict
                 ELSE
                    acceff=zero
                 ENDIF
                 pos(p,k)=pos(p,k)+vel(p,k)*tsteppos*(one+acceff)
 500          CONTINUE
 600       CONTINUE

        ENDIF

C   Update position time, system time.
C   ----------------------------------
        tpos=tpos+tsteppos
        tnow=tpos

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE vextrap
C
C
C***********************************************************************
C
C
C     Subroutine to extrapolate velocities forward in time using 
C     current estimate of acceleration.  Primarily used to obtain
C     an estimate of the velocity which is time-synchronized with
C     the positions, for the SPH particles.  The values of veltpos 
C     are then used in computing velocity-dependent forces, to
C     maintain time-centering when updating vel.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k

C=======================================================================

C   Loop over all velocity components for all SPH particles.
C   --------------------------------------------------------


           DO 200 k=1,ndim
              DO 100 p=1,nsph
                 veltpos(p,k)=vel(p,k)+acc(p,k)*(tnow-tvel(p))
 100          CONTINUE
 200       CONTINUE


        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE stepvel
C
C
C***********************************************************************
C
C
C     Subroutine to advance the velocities of the bodies for a
C     timestep dt.  This version allows for an optional frictional 
C     damping term in the equations of motion for the SPH particles.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k,i
        REAL acceff,vfac1now,vfac2now,dt

C=======================================================================

C   Loop over all velocity components for all bodies.
C   -------------------------------------------------

        IF(.NOT.friction) THEN

           IF(.NOT.comove) THEN

              DO 200 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 100 i=1,npactive
                    p=pactive(i)
                    vel(p,k)=vel(p,k)+acc(p,k)*dtime/itimestp(p)
 100             CONTINUE

 200          CONTINUE

           ELSE

              DO 250 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP

                 DO 220 i=1,npactive
                    p=pactive(i)
                    dt=dtime/itimestp(p)
                    vfac1now=(one-hubble*dt)/(one+hubble*dt)
                    vfac2now=dt/(cosmof3*(one+hubble*dt))
                    vel(p,k)=vel(p,k)*vfac1now+acc(p,k)*vfac2now
 220             CONTINUE

 250          CONTINUE

           ENDIF

        ELSE

           IF(.NOT.friclucy) THEN

              DO 400 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 300 i=1,npactive
                    p=pactive(i)
                    dt=dtime/itimestp(p)
                    IF(p.LE.nsph) THEN
                       acceff=cfrict/dt
                    ELSE
                       acceff=zero
                    ENDIF 
                    vel(p,k)=(vel(p,k)*(one-0.5*dt*acceff)+acc(p,k)*dt)/
     &                       (one+0.5*dt*acceff)
 300             CONTINUE

 400          CONTINUE

           ELSE

              DO 600 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 500 i=1,npactive
                    p=pactive(i)
                    dt=dtime/itimestp(p)
                    IF(p.LE.nsph) THEN
                       acceff=cfrict
                    ELSE
                       acceff=zero
                    ENDIF 
                    vel(p,k)=(vel(p,k)+acc(p,k)*dt)/(one+acceff)
 500             CONTINUE

 600          CONTINUE

           ENDIF

        ENDIF

C   Update velocity times.
C   ----------------------
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 700 i=1,nsphact
           p=pactive(i)
           tvel(p)=tvel(p)+dtime/itimestp(p)
 700    CONTINUE        

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE vfilter
C
C
C***********************************************************************
C
C
C     Subroutine to filter velocities of SPH particles according to:
C
C             v-filter = (1 - epsfiltr) * v + epsfiltr * <v>
C
C     where <v> is a smoothed estimate of the velocity field:
C                    
C             <v>_i = SUM_j { m_j v_j W_ij} / SUM_j {m_j W_ij}  (1)
C
C     The kernel used in (1) is the spherically symmetric M-5 kernel,
C     as opposed to the M-4 kernel throughout the rest of the code.
C     The M-5 kernel is, in 3-D:
C
C               /          
C               | 115/8 - 15v**2 + 6 v**4                0 <= v <= 1/2
C               |              
C               | 55/4 + 5v - 30v**2 + 20 v**3 - 4v**4   1/2 <= v <= 3/2
C   W_5 = wnorm |
C               | (5/2 - v)**4                           3/2 <= v <= 5/2
C               |  
C               |  0                                      v => 5/2
C               \
C
C     The factor wnorm is a normalization constant and is equal to
C     wnorm = 1./(20 * pi * h**3).  For present purposes, the kernel
C     must be rescaled to the interval 0 <= v <= 2.  This gives:
C
C               /          
C               | 115/8 - 15*25v**2/16 + 6*625v**4/256   0 <= v <= 2/5
C               |              
C               | 55/4 + 25v/4 - 30*25v**2/16 + 
C               |        20*125*v**3/64 - 625v**4/64     2/5 <= v <= 6/5
C W_5s = wnorms |
C               | (5/2 - 5*v/4)**4                       6/5 <= v <= 2
C               |  
C               |  0                                      v => 2
C               \
C
C
C     where wnorms = 25/(256 * pi * h**3).
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),nnind(nbodsmax),inear,iwsm,iwsm1,
     &          jsubset(nbodsmax)
        LOGICAL firstc
        REAL dx,dy,dz,dr2,hsminv,drw,distnorm,dr2p,dr2i,drw1,
     &       w5smooth(0:1+ninterp),deldr2,xw,xw2,vnorm(nworkvec),
     &       wmassavg,wnorm,wsm,wmass,wsm1,wmass1,vnormp,sphvxp,
     &       sphvyp,sphvzp

        EQUIVALENCE (nb(1),bodlist(1)),(nnind(1),templist(1)),
     &              (jsubset(1),subindex(1)),(vnorm(1),workvect(1))

        DATA firstc /.TRUE./

        SAVE firstc,w5smooth

C=======================================================================

C   Initialize properties of smoothing kernel.
C   ------------------------------------------

        IF(firstc) THEN

           firstc=.FALSE.

           deldr2=4./ninterp

           DO 10 i=0,1+ninterp
              xw2=i*deldr2
              xw=SQRT(xw2)

              IF(xw2.LE.0.4*0.4) THEN
                 w5smooth(i)=25.*(115./8.-15.*25.*xw2/16.+
     &                       6.*625.*xw2*xw2/256.)/256.            
              ENDIF

              IF(xw2.GT.0.4*0.4.AND.xw2.LE.1.2*1.2) THEN
                 w5smooth(i)=25.*(55./4.+25.*xw/4.-30.*25.*xw2/16.+
     &                       20.*125.*xw*xw2/64.-625.*xw2*xw2/64.)/256.
              ENDIF

              IF(xw2.GT.1.2*1.2.AND.xw2.LE.four) THEN
                 w5smooth(i)=25.*(2.5-1.25*xw)**4/256.
              ENDIF

              IF(xw2.GT.four) w5smooth(i)=zero

 10        CONTINUE

        ENDIF


        DO 20 p=1,nsph
           isubset(p)=0
           vnorm(p)=0.
           tempsphv(p,1)=0.
           tempsphv(p,2)=0.
           tempsphv(p,3)=0.
 20     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 p=1,nsphact
           isubset(pactive(p))=1
 30     CONTINUE

C   Compute filter properties, maintaining inversion symmetry.
C   ----------------------------------------------------------
        DO 90 p=1,nsph

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 35 i=1,nnearlis(p)
              jsubset(i)=isubset(nearbods(pnear(p)+i-1))+isubset(p)
 35        CONTINUE

           CALL WHENNE(nnearlis(p),jsubset,1,0,nnind,inear)

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

           vnormp=0.0
           sphvxp=0.0
           sphvyp=0.0
           sphvzp=0.0

           IF(symmetry.EQ.'hk') THEN

              hsminv=1.0/hsmooth(p)
              wnorm=piinv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 40 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 wsm=(1.-drw)*w5smooth(iwsm)+drw*w5smooth(1+iwsm)
                 wmass=wnorm*wsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 wsm1=(1.-drw1)*w5smooth(iwsm1)+drw1*w5smooth(1+iwsm1)
                 wmass1=piinv*wsm1/(hsmooth(nb(i))*hsmooth(nb(i))**2)

                 wmassavg=0.5*(wmass+wmass1)*isubset(p)

                 vnormp=vnormp+mass(nb(i))*wmassavg
                 sphvxp=sphvxp+mass(nb(i))*vel(nb(i),1)*wmassavg
                 sphvyp=sphvyp+mass(nb(i))*vel(nb(i),2)*wmassavg
                 sphvzp=sphvzp+mass(nb(i))*vel(nb(i),3)*wmassavg

                 wmassavg=0.5*(wmass+wmass1)*isubset(nb(i))

                 vnorm(nb(i))=vnorm(nb(i))+mass(p)*wmassavg
                 tempsphv(nb(i),1)=tempsphv(nb(i),1)+mass(p)*
     &                             vel(p,1)*wmassavg
                 tempsphv(nb(i),2)=tempsphv(nb(i),2)+mass(p)*
     &                              vel(p,2)*wmassavg
                 tempsphv(nb(i),3)=tempsphv(nb(i),3)+mass(p)*
     &                             vel(p,3)*wmassavg

 40           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 50 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 wnorm=piinv*hsminv**2*hsminv**2
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 wsm=(1.-drw)*w5smooth(iwsm)+drw*w5smooth(1+iwsm)

                 wmass=wnorm*wsm*isubset(p)

                 vnormp=vnormp+mass(nb(i))*wmass
                 sphvxp=sphvxp+mass(nb(i))*vel(nb(i),1)*wmass
                 sphvyp=sphvyp+mass(nb(i))*vel(nb(i),2)*wmass
                 sphvzp=sphvzp+mass(nb(i))*vel(nb(i),3)*wmass

                 wmass=wnorm*wsm*isubset(nb(i))

                 vnorm(nb(i))=vnorm(nb(i))+mass(p)*wmass
                 tempsphv(nb(i),1)=tempsphv(nb(i),1)+mass(p)*
     &                             vel(p,1)*wmass
                 tempsphv(nb(i),2)=tempsphv(nb(i),2)+mass(p)*
     &                             vel(p,2)*wmass
                 tempsphv(nb(i),3)=tempsphv(nb(i),3)+mass(p)*
     &                             vel(p,3)*wmass

 50           CONTINUE

           ENDIF

C   Add self-interaction term.
C   --------------------------

           IF(isubset(p).EQ.1) THEN

              IF(symmetry.EQ.'be') THEN
                 wnorm=piinv/hsmooth(p)**3
              ENDIF

              wmass=wnorm*mass(p)*25.*115./(256.*8.)
              vnorm(p)=vnorm(p)+wmass+vnormp
              tempsphv(p,1)=tempsphv(p,1)+wmass*vel(p,1)+sphvxp
              tempsphv(p,2)=tempsphv(p,2)+wmass*vel(p,2)+sphvyp
              tempsphv(p,3)=tempsphv(p,3)+wmass*vel(p,3)+sphvzp

           ENDIF

 90     CONTINUE

C   Apply filter to SPH velocities.
C   -------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 100 i=1,nsphact
           p=pactive(i)
           vel(p,1)=(one-epsfiltr)*vel(p,1)+epsfiltr*tempsphv(p,1)/
     &              vnorm(p)
           vel(p,2)=(one-epsfiltr)*vel(p,2)+epsfiltr*tempsphv(p,2)/
     &              vnorm(p)
           vel(p,3)=(one-epsfiltr)*vel(p,3)+epsfiltr*tempsphv(p,3)/
     &              vnorm(p)
 100    CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE stepsph
C
C
C***********************************************************************
C
C
C     Subroutine to advance the thermal energy of the SPH particles
C     and compute the hydrodynamic contribution to the acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C=======================================================================

        CALL maketsph
C            --------
        CALL stepnear
C            --------
        CALL density
C            -------

        IF(artfvisc.EQ.'bulk') THEN
 
           IF(.NOT.isotherm) THEN

              IF(.NOT.uentropy) THEN

                 CALL stepeth('predict')
C                     --------
                 CALL ethdotbv('predict')
C                     --------
                 CALL stepeth('correct')
C                     --------
                 CALL ethdotbv('correct')
C                     --------
              ELSE

                 CALL stepent('predict')
C                     --------
                 CALL entdotbv('predict')
C                     --------
                 CALL stepent('correct')
C                     --------
                 CALL entdotbv('correct')
C                     --------
              ENDIF
           ENDIF

           IF(.NOT.endstep) CALL accsphbv
C                                --------
        ELSE
           IF(artfvisc.EQ.'sphv') THEN

              IF(.NOT.isotherm) THEN

                 IF(.NOT.uentropy) THEN
    
                    CALL stepeth('predict')
C                        --------
                    CALL ethdotcv('predict')
C                        --------
                    CALL stepeth('correct')
C                        --------
                    CALL ethdotcv('correct')
C                        --------
                 ELSE
                    CALL stepent('predict')
C                        --------
                    CALL entdotcv('predict')
C                        --------
                    CALL stepent('correct')
C                        --------
                    CALL entdotcv('correct')
C                        --------
                 ENDIF
              ENDIF

              IF(.NOT.endstep) CALL accsphcv
C                                   --------
           ELSE
              IF(.NOT.isotherm) THEN
   
                 IF(.NOT.uentropy) THEN

                    CALL stepeth('predict')
C                        --------
                    CALL ethdot('predict')
C                        ------
                    CALL stepeth('correct')
C                        --------
                    CALL ethdot('correct')
C                        ------
                 ELSE

                    CALL stepent('predict')
C                        --------
                    CALL entdot('predict')
C                        ------
                    CALL stepent('correct')
C                        --------
                    CALL entdot('correct')
C                        ------
                 ENDIF
              ENDIF

              IF(.NOT.endstep) CALL accsph
C                                   ------
           ENDIF

        ENDIF

        IF(isotherm) CALL isomumax
C                         --------

        IF(radiate.AND.(.NOT.isotherm)) CALL ethrad
C                                            ------

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stepnear
C
C
C***********************************************************************
C
C
C     Subroutine to update smoothing lengths and establish nearest 
C     neighbor lists for all SPH particles.  The list of neighbors
C     is returned from the subroutine findnear in the vector 
C     nearlist, which is equivalenced to the common array bodlist.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,pbody,npnear,inear,insave,nnearbod,nearpb(nbodsmax),
     &          np,jsubset(nbodsmax),nearlist(nbodsmax),
     &          keepnear(nbodsmax)
        LOGICAL keepcrit,savecrit
        REAL testnn,fourhsm2,dx,dy,dz

        EQUIVALENCE (nearlist(1),bodlist(1)),(jsubset(1),subindex(1)),
     &              (nearpb(1),parent(1)),(keepnear(1),asubp(1))

C=======================================================================

C   If required, extrapolate smoothing lengths.
C   -------------------------------------------
        IF(variablh) CALL hextrap
C                         -------

C   Initialize nearest neighbor diagnostics.
C   ----------------------------------------
        nnearbod=0
        nntot=0
        nnmin=nsphmax
        nnmax=0

C   Find nearest neighbors of grouped cells.
C   ----------------------------------------

        DO 100 p=1,ngroups

C   Find nearest neighbors.
C   -----------------------
           CALL findnear(p,npnear)
C               --------

           DO 70 np=pgroupb(p),pgroupb(p+1)-1

              pbody=groupbod(np)

              fourhsm2=four*hsmooth(pbody)*hsmooth(pbody)

              DO 20 i=1,npnear
                 dx=pos(pbody,1)-pos(nearlist(i),1)
                 dy=pos(pbody,2)-pos(nearlist(i),2)
                 dz=pos(pbody,3)-pos(nearlist(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 tempvect(i)=dx**2+dy**2+dz**2
                 IF(nearlist(i).EQ.pbody) tempvect(i)=fourhsm2+one
 20           CONTINUE

C   Filter out neighbors further than two smoothing lengths.
C   --------------------------------------------------------

              CALL WHENFLT(npnear,tempvect,1,fourhsm2,isubset,inear)

              testnn=ABS(REAL(inear-nsmooth)/REAL(nsmooth))
              savecrit=(.NOT.variablh).OR.(testnn.LE.nsmtol).OR.
     &                 (templist(pbody).GT.2.AND.inear.GT.nsmooth)

              IF(savecrit) THEN

                 templist(pbody)=-1

                 nntot=nntot+inear
                 IF(inear.LT.nnmin) nnmin=inear
                 IF(inear.GT.nnmax) nnmax=inear
                 nnear(pbody)=inear

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 30 i=1,inear
                    jsubset(i)=nearlist(isubset(i))
                    keepcrit=hsmooth(pbody).GT.hsmooth(jsubset(i)).OR.
     &                       ((hsmooth(pbody).EQ.hsmooth(jsubset(i)))
     &                       .AND.pbody.LT.jsubset(i))
                    IF(keepcrit) THEN
                       keepnear(i)=2
                    ELSE
                       keepnear(i)= -2
                    ENDIF
 30              CONTINUE
 
                 CALL WHENIGT(inear,keepnear,1,0,jsubset,insave)

                 pnear(pbody)=nnearbod+1
                 nnearlis(pbody)=insave
                 nnearbod=nnearbod+insave

                 IF(nnearbod.GT.maxnearb)
     &              CALL terror(' array overflow in stepnear ')
C                        ------

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 40 i=1,insave
                    nearbods(nnearbod-insave+i)=
     &                       nearlist(isubset(jsubset(i)))
 40              CONTINUE

              ELSE

                 IF(inear.GT.nsmooth) THEN

                    DO 50 i=1,inear
                       nearpb(i)=nearlist(isubset(i))
 50                 CONTINUE
 
                    CALL reduceh(pbody,inear,nnearbod)
C                        -------

                 ELSE

                    CALL enlargeh(pbody,nnearbod)
C                        --------
                 ENDIF

              ENDIF

 70        CONTINUE

 100    CONTINUE

        nnavg=nntot/nsph

        IF(variablg) THEN

           DO 200 i=1,nsph
              epsvect(i)=hsmooth(i)
 200       CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE hextrap
C
C
C***********************************************************************
C
C
C     Subroutine to extrapolate smoothing lengths forward in time
C     by requiring that each particle interact with roughly a 
C     constant number of near neighbors.  The extrapolation is
C     performed in accord with the Courant condition.  The variable
C     templist is used to return the outcome of the Courant test.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,nnearfix
        REAL hsmc,dt,hsmctemp,abshdivv

C=======================================================================

        DO 10 p=1,nsph

           dt=tnow-teth

           IF(nnear(p).EQ.0) THEN
              nnearfix=1
           ELSE
              nnearfix=0
           ENDIF
           hsmooth(p)=hsmooth(p)*0.5*((REAL(nsmooth)/REAL(nnear(p)+
     &                nnearfix))**(one/3.)+one)
           abshdivv=hsmdivv(p)
           abshdivv=ABS(abshdivv)
           hsmc=dt*(csound(p)+1.2*(alpha*csound(p)+beta*mumaxdvh(p))+
     &          abshdivv)/courant
           hsmc=hsmc/cosmofac
           IF(hsmc.LT.hsmooth(p)) THEN
              hsmctemp=one
           ELSE
              hsmctemp=four
           ENDIF
           templist(p)=INT(hsmctemp)
           IF(hsmooth(p).LT.hsmc) hsmooth(p)=hsmc
 10     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE density
C
C
C***********************************************************************
C
C
C     Subroutine to compute smoothed estimate of the local density
C     and divergence and magnitude of the curl of the velocity field 
C     for SPH particles.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),iwsm,iwsm1
        REAL hsminv,wnorm,distnorm,dr2,dr2p,drw,wsm,wmass,dr2i,drw1,
     &       wsm1,wmass1,dwnorm,vdotdr,dwsm,dwmass,dwsm1,dwmass1,dx,dy,
     &       dz,dvx,dvy,dvz,curlvx(nsphmax),curlvy(nsphmax),
     &       curlvz(nsphmax),dwmp,dwmnbi

        EQUIVALENCE (nb(1),bodlist(1)),(curlvx(1),tempsphv(1,1)),
     &              (curlvy(1),tempsphv(1,2)),(curlvz(1),tempsphv(1,3))

C=======================================================================

C   Zero the density and velocity divergence.
C   -----------------------------------------
        DO 10 p=1,nsph
           rho(p)=0.
           hsmdivv(p)=0.
           curlvx(p)=0.
           curlvy(p)=0.
           curlvz(p)=0.
 10     CONTINUE

C   Compute density and div v of particles.
C   ---------------------------------------

        DO 30 p=1,nsph


              hsminv=1.0/hsmooth(p)
              wnorm=piinv*hsminv*hsminv*hsminv
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 20 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                 wmass=wnorm*wsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 wsm1=(1.-drw1)*wsmooth(iwsm1)+drw1*wsmooth(1+iwsm1)
                 wmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                  hsmooth(nb(i)))*wsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 dvx=veltpos(p,1)-veltpos(nb(i),1)+dx*cosmohub
                 dvy=veltpos(p,2)-veltpos(nb(i),2)+dy*cosmohub
                 dvz=veltpos(p,3)-veltpos(nb(i),3)+dz*cosmohub
                 vdotdr=dvx*dx+dvy*dy+dvz*dz
                 dwmp=mass(p)*0.5*(dwmass+dwmass1)
                 dwmnbi=mass(nb(i))*0.5*(dwmass+dwmass1)
                 rho(p)=rho(p)+0.5*mass(nb(i))*(wmass+wmass1)
                 hsmdivv(p)=hsmdivv(p)-dwmnbi*vdotdr
                 curlvx(p)=curlvx(p)+dwmnbi*(dz*dvy-dy*dvz)
                 curlvy(p)=curlvy(p)+dwmnbi*(dx*dvz-dz*dvx)
                 curlvz(p)=curlvz(p)+dwmnbi*(dy*dvx-dx*dvy)
                 rho(nb(i))=rho(nb(i))+0.5*mass(p)*(wmass+wmass1)
                 hsmdivv(nb(i))=hsmdivv(nb(i))-dwmp*vdotdr
                 curlvx(nb(i))=curlvx(nb(i))+dwmp*(dz*dvy-dy*dvz)
                 curlvy(nb(i))=curlvy(nb(i))+dwmp*(dx*dvz-dz*dvx)
                 curlvz(nb(i))=curlvz(nb(i))+dwmp*(dy*dvx-dx*dvy)
 20           CONTINUE

           wmass=wnorm*mass(p)
           rho(p)=rho(p)+wmass

 30     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 p=1,nsph
           hsmdivv(p)=hsmooth(p)*hsmdivv(p)/rho(p)
           hsmcurlv(p)=hsmooth(p)*SQRT(curlvx(p)**2+curlvy(p)**2+
     &                 curlvz(p)**2)/rho(p)
 40     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE stepeth(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the thermal energy implicitly by a
C     timestep dt, using the Newton-Raphson method.  The local
C     variable p is a pointer to the bodies.
C
C     In this version of stepeth, and in the subroutines ethrad and
C     ethdiff (the latter called by bisect), radiative cooling is cut
C     off below a critical temperature, tcutoff, determined by a 
C     Jeans mass criterion.  Specifically, if the Jeans mass is 
C     written in the form:
C
C        M-J = mjconst * csound**4 / SQRT(P * G**3)
C
C     then for a given M-J the cutoff temperature is given by:
C
C        tcutoff = ctcutoff * meanmwt * mhboltz * G * 
C                  (M-J)**(2/3) * rho **(1/3)
C
C     where meanmwt and mhboltz are as defined in the main routine,
C     G = 1 in the adopted system of units.  For the classical 
C     definition of the Jeans mass, mjconst = pi **(3/2) while for 
C     the Bonnor-Ebert Jeans mass, mjconst = 1.18.  The constant 
C     ctcutoff, input via the parameter file, includes not only 
C     (mjconst)**(-2/3), but a factor determining the number of 
C     particles to include in the Jeans mass.  For the two choices 
C     above, (mjconst)**(-2/3) = 0.31831 and 0.89553, respectively.
C     Typically, the factor from the Jeans mass will be ~ 1 --> 
C     (nsmooth/8) ** (2/3), implying a total variation in ctcutoff 
C     of ~ 10.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        CHARACTER*7 pc
        INTEGER p,nit,i,ib,nitmax,nsubset,nblist,ifirst
        REAL dt,smallnum,temp,templog,tenf1,tenf2,heatcool,dheatc,
     &       tempmin,tempmax,ethlast,ethmin,ethmax,ethtemp,tempdiff,
     &       cutoff,tenf3,tenf3exp,compcool,aconst,tcutoff,
     &       tcutoffc,compdiff,cutoffc
     
        PARAMETER(nitmax=15)
     
        SAVE smallnum,tempmin,tempmax,ifirst,ethmin,ethmax
     
        DATA smallnum/1.e-10/,tempmin/10.0/,tempmax/1.e10/,ifirst/1/
     
C=======================================================================
     
        dt=tnow-teth
     
        IF(ifirst.EQ.1) THEN
           ethmin=tempmin/(gamma1*meanmwt*mhboltz)
           ethmax=tempmax/(gamma1*meanmwt*mhboltz)
           ifirst=0
        ENDIF
     
        IF(pc.EQ.'predict') THEN
           DO 5 p=1,nsph
              dethold(p)=dethdt(p)
              dethdt(p)=cosmfold*dethdt(p)/cosmofac
 5         CONTINUE
        ENDIF

        DO 10 p=1,nsph
           ethold(p)=ethermal(p)
           bodlist(p)=p
 10     CONTINUE

        nit=0
        nblist=nsph
     
 20     CONTINUE
     
        IF(.NOT.radiate) THEN
     
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nblist
              ib=bodlist(i)
              ethermal(ib)=ethold(ib)+0.5*dt*(dethdt(ib)+dethold(ib))
              ethtemp=ethermal(ib)
              IF(ethtemp.LE.zero) THEN
                 tempvect(i)=two*smallnum
              ELSE
                 tempvect(i)=zero
              ENDIF
 30        CONTINUE
     
        ELSE
     
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nblist
              ib=bodlist(i)
              temp=meanmwt*mhboltz*gamma1*ethermal(ib)
              tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(ib)*nsmooth/8./
     &                (rho(ib)/cosmof3))**(2./3.))-((epsvect(ib)/
     &                piinv)**2.)/3.)*(rho(ib)/cosmof3)
              IF(tcutoff.LT.mintemp) tcutoff=mintemp
              templog=LOG10(temp)
              tempdiff= 30.*(1.-temp/tcutoff)
              IF(tempdiff.LT.0.)tempdiff=0.
              cutoff=exp(-tempdiff *tempdiff)

              tcutoffc=tcutoff
              IF(tcutoff.LT.10000.) tcutoffc=10000.
              compdiff= 30.*(1.-temp/tcutoffc)
              IF(compdiff.LE.0.)compdiff=0.
              cutoffc=exp(-compdiff *compdiff)

              tenf1=10.**(-4.43-0.273*(4.-templog)**2)
              tenf2=10.**(-0.10-1.880*(5.23-templog)**4)
              tenf3exp=6.2-templog
              IF(templog-6.2.GT.zero) tenf3exp=zero
              tenf3=10.**(-2.*tenf3exp**4-1.7)
              compcool=comptmh*((1.+redshift)**4)

              heatcool=aheatmh*fhydrogn+bheatmh2*rho(ib)*fhydrog2/
     &                 cosmof3-ccoolmh2*rho(ib)*(tenf1+tenf2+tenf3)*
     &                 cutoff*fhydrog2/cosmof3-compcool*temp*cutoffc

              aconst=slowcool*ethermal(ib)/dt+dethdt(ib)
              IF(aconst.LE.zero) aconst=tiny
              aconst=heatcool/aconst
              heatcool=heatcool/SQRT(1.+aconst*aconst)
              tempvect(i)=ethermal(ib)-ethold(ib)-0.5*dt*(dethdt(ib)+
     &                    heatcool+dethold(ib)+derad(ib))

              dheatc= - ccoolmh2*rho(ib)/cosmof3*(((0.546*(4.-
     &                templog)*tenf1+
     &                7.52*((5.23-templog)**3)*tenf2+8.*tenf3exp**3*
     &                tenf3)*cutoff/ethermal(ib))+((tenf1+tenf2+tenf3)*
     &                2.*30.*cutoff*tempdiff*
     &                meanmwt*mhboltz*gamma1/tcutoff))*fhydrog2-
     &                compcool*meanmwt*mhboltz*gamma1*
     &                (cutoffc+compcool*temp*2.*30.*cutoffc*compdiff/
     &                tcutoffc)

              dheatc=(dheatc*(1.-aconst*aconst/(1.+aconst*aconst))+
     &                aconst*aconst*aconst/(1.+aconst*aconst)*slowcool
     &                /dt)/SQRT(1.+aconst*aconst)
              tempvect(i)=tempvect(i)/(1.-0.5*dt*dheatc)
 40        CONTINUE
     
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nblist
              ib=bodlist(i)
              ethlast=ethermal(ib)
              ethermal(ib)=ethermal(ib)-tempvect(i)
              ethtemp=ethermal(ib)
              IF(ethmax.LT.ethtemp) ethtemp=ethmax
              IF(ethtemp.LT.ethmin) ethtemp=ethmin
              ethermal(ib)=ethtemp
              tempvect(i)=(ethlast-ethermal(ib))/ethlast
              tempvect(i)=ABS(tempvect(i))
              IF(ethtemp.LT.1.05*ethmin.OR.ethtemp.LT.ethold(ib)/10.) 
     &                  tempvect(i)=two*smallnum
              IF(ethtemp.GT.0.95*ethmax.OR.ethtemp.GT.10.*ethold(ib)) 
     &                  tempvect(i)=two*smallnum
 50        CONTINUE
     
        ENDIF
     
        nit = nit + 1
     
        CALL WHENFGT(nblist,tempvect,1,smallnum,isubset,nsubset)
     
        DO 60 i=1,nsubset
           templist(i)=bodlist(isubset(i))
 60     CONTINUE
     
        DO 70 i=1,nsubset
           bodlist(i)=templist(i)
 70     CONTINUE
     
        nblist=nsubset
     
        IF(nit.LE.nitmax) THEN
     
           IF(nblist.GT.0) GO TO 20
     
        ELSE

           CALL bisect(nblist,tempmin,tempmax,smallnum,dt)
C               ------     
     
        ENDIF
     
        IF(radiate) CALL ethcheck(dt,smallnum,tempmin,tempmax)
C                           --------
C   Recompute sound speed, pressure / rho**2.
C   -----------------------------------------
        DO 80 p=1,nsph
           csound(p)=SQRT(gamma*gamma1*ethermal(p))
 80     CONTINUE
     
        IF(pc.EQ.'correct') THEN
           teth=teth+dt
        ENDIF

        IF(pc.EQ.'predict') THEN
           DO 90 p=1,nsph
              ethermal(p)=ethold(p)
 90        CONTINUE
        ENDIF
     
        RETURN
        END
C***********************************************************************
C
C
                     SUBROUTINE ethcheck(dt,smallnum,tempmin,tempmax)
C
C
C***********************************************************************
C
C
C     Subroutine to check whether the new value for the thermal 
C     energy is actually a root.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER ib,nerror
        REAL dt,tenf1,tenf2,heatcool,tempdiff,cutoff,tenf3,
     &       tenf3exp,temp,logtemp,compcool,aconst,tcutoff,
     &       smallnum,tcutoffc,compdiff,cutoffc,tempmin,tempmax

C=======================================================================
        DO 10 ib=1,nsph

           IF(ethermal(ib)*meanmwt*mhboltz*gamma1.LT.1.05*tempmin)THEN
              tempvect(ib)=smallnum/10.
           ELSE
              temp=ethermal(ib)*meanmwt*mhboltz*gamma1
              logtemp=LOG10(temp)

              tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(ib)*nsmooth/8./
     &                (rho(ib)/cosmof3))**(2./3.))-((epsvect(ib)/
     &                piinv)**2.)/3.)*(rho(ib)/cosmof3)
              IF(tcutoff.LT.mintemp) tcutoff=mintemp
              tempdiff= 30.*(1.-temp/tcutoff)
              IF(tempdiff.LE.0.)tempdiff=0.
              cutoff=exp(-tempdiff *tempdiff)

              tcutoffc=tcutoff
              IF(tcutoff.LT.10000.) tcutoffc=10000.
              compdiff= 30.*(1.-temp/tcutoffc)
              IF(compdiff.LE.0.)compdiff=0.
              cutoffc=exp(-compdiff *compdiff)

              tenf1=10.**(-4.43-0.273*(4.-logtemp)**2)
              tenf2=10.**(-0.10-1.880*(5.23-logtemp)**4)
              tenf3exp=6.2-logtemp
              IF(logtemp-6.2.GT.zero) tenf3exp=zero
              tenf3=10.**(-2.*tenf3exp**4-1.7)
              compcool=comptmh*temp*((1.+redshift)**4)

              heatcool=aheatmh*fhydrogn+bheatmh2*rho(ib)*fhydrog2/
     &                 cosmof3-ccoolmh2*rho(ib)*(tenf1+tenf2+tenf3)*
     &                 cutoff*fhydrog2/cosmof3-compcool*cutoffc

              aconst=slowcool*temp/(meanmwt*mhboltz*gamma1*dt)+
     &               dethdt(ib)
              IF(aconst.LE.zero) aconst=tiny
              aconst=heatcool/aconst
              heatcool=heatcool/SQRT(1.+aconst*aconst)
              tempvect(ib)=abs(temp/(meanmwt*mhboltz*gamma1)-
     &                ethold(ib)-0.5*dt*
     &                (dethdt(ib)+heatcool+dethold(ib)+derad(ib)))
           ENDIF

 10     CONTINUE

        CALL WHENFGT(nsph,tempvect,1,10.*smallnum,bodlist,nerror)

        IF(nerror.GT.0)THEN
           CALL bisect(nerror,tempmin,tempmax,smallnum,dt)
C               ------
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE ethdot(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the time derivative of the thermal
C     energy.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,nb(nbodsmax),iwsm,imumax,ismax,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,muij(nbodsmax),dwsm,vdotdr,
     &       piij,tmuij,dr2p,distnorm,dwmass,dr2i,drw1,dwsm1,dwmass1,
     &       hsmavg,eijp,eijnbi,rhoprhoi

        EQUIVALENCE (nb(1),bodlist(1)),(muij(1),tempvect(1))

C=======================================================================

C   Zero out derivative of thermal energy.
C   --------------------------------------
        DO 10 p=1,nsph
           dethdt(p)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute derivative of thermal energy.
C-----------------------------------------------------------------------

        DO 40 p=1,nsph

           IF(symmetry.EQ.'hk') THEN
              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    eijp=vdotdr*(csound(nb(i))*csound(p)/
     &                   (SQRT(rhoprhoi)*gamma)+0.5*piij)
                    eijnbi=eijp
                 ELSE
                    eijp=vdotdr*(csound(p)**2/(gamma*rho(p))+0.5*piij)
                    eijnbi=vdotdr*(csound(nb(i))**2/(gamma*rho(nb(i)))+
     &                   0.5*piij)
                 ENDIF
                 dethdt(p)=dethdt(p)+eijp*mass(nb(i))*0.5*(dwmass+
     &                     dwmass1)
                 dethdt(nb(i))=dethdt(nb(i))+eijnbi*mass(p)*0.5*(dwmass+
     &                         dwmass1)
                 muij(i)=ABS(muij(i))
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    eijp=vdotdr*(csound(nb(i))*csound(p)/
     &                   (SQRT(rhoprhoi)*gamma)+0.5*piij)
                    eijnbi=eijp
                 ELSE
                    eijp=vdotdr*(csound(p)**2/(gamma*rho(p))+0.5*piij)
                    eijnbi=vdotdr*(csound(nb(i))**2/(gamma*rho(nb(i)))+
     &                   0.5*piij)
                 ENDIF
                 dethdt(p)=dethdt(p)+eijp*mass(nb(i))*dwmass
                 dethdt(nb(i))=dethdt(nb(i))+eijnbi*mass(p)*dwmass
                 muij(i)=ABS(muij(i))
 35           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0.AND.pc.EQ.'correct') THEN
              imumax=ISMAX(nnearlis(p),muij,1)
              tmuij=muij(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 60 i=1,nnearlis(p)
                 tmuij=muij(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 60           CONTINUE

           ENDIF

 40     CONTINUE

C##### LH 050990
 
        DO 70 p=1,nsph
           dethdt(p)=dethdt(p)/cosmofac
 70     CONTINUE
		      
C##### LH 050990

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE entdot(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the time derivative of the entropic
C     function a(s).
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,nb(nbodsmax),iwsm,imumax,ismax,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,muij(nbodsmax),dwsm,vdotdr,
     &       piij,tmuij,dr2p,distnorm,dwmass,dr2i,drw1,dwsm1,dwmass1,
     &       hsmavg,eij

        EQUIVALENCE (nb(1),bodlist(1)),(muij(1),tempvect(1))

C=======================================================================

C   Zero out derivative of thermal energy.
C   --------------------------------------
        DO 10 p=1,nsph
           dentdt(p)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute derivative of thermal energy.
C-----------------------------------------------------------------------

        DO 40 p=1,nsph

           IF(symmetry.EQ.'hk') THEN

              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 eij=rho(p)*rho(nb(i))
                 eij=vdotdr*0.5*piij
                 dentdt(p)=dentdt(p)+eij*mass(nb(i))*0.5*(dwmass+
     &                     dwmass1)
                 dentdt(nb(i))=dentdt(nb(i))+eij*mass(p)*0.5*(dwmass+
     &                         dwmass1)
                 muij(i)=ABS(muij(i))
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 eij=rho(p)*rho(nb(i))
                 eij=vdotdr*0.5*piij
                 dentdt(p)=dentdt(p)+eij*mass(nb(i))*dwmass
                 dentdt(nb(i))=dentdt(nb(i))+eij*mass(p)*dwmass
                 muij(i)=ABS(muij(i))
 35           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0.AND.pc.EQ.'correct') THEN
              imumax=ISMAX(nnearlis(p),muij,1)
              tmuij=muij(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 60 i=1,nnearlis(p)
                 tmuij=muij(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 60           CONTINUE

           ENDIF

 40     CONTINUE

C##### LH 050990

        DO 70 p=1,nsph
           dentdt(p)=dentdt(p)*gamma1/(cosmofac*rho(p)**gamma1)
 70     CONTINUE

C##### LH 050990

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE accsph
C
C
C***********************************************************************
C
C
C     Subroutine to compute SPH contribution to the acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),nnind(nbodsmax),inear,iwsm,
     &          j,jsubset(nbodsmax),ninear,iwsm1
        REAL dx,dy,dz,hsmavg,dr2,hsminv,drw,dwnorm,dwsm1,muij(nbodsmax),
     &       dwsm,vdotdr,aij(nworkvec),piij,distnorm,dr2p,dwmass,
     &       dr2i,drw1,dwmass1,aijmass,rhoprhoi,gradprho

        EQUIVALENCE (nb(1),bodlist(1)),(nnind(1),templist(1)),
     &              (jsubset(1),subindex(1)),(muij(1),tempvect(1)),
     &              (aij(1),workvect(1))

C=======================================================================

        DO 10 p=1,nsph
           isubset(p)=0
 10     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 20 p=1,nsphact
           isubset(pactive(p))=1
 20     CONTINUE

C   Compute SPH acceleration, maintaining inversion symmetry.
C   ---------------------------------------------------------
        DO 90 p=1,nsph

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 25 i=1,nnearlis(p)
              jsubset(i)=isubset(nearbods(pnear(p)+i-1))+isubset(p)
 25        CONTINUE

           CALL WHENNE(nnearlis(p),jsubset,1,0,nnind,inear)

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

           IF(symmetry.EQ.'hk') THEN

              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    gradprho=2.*csound(p)*csound(nb(i))/(gamma*
     &                       SQRT(rhoprhoi))
                 ELSE
                    gradprho=csound(p)**2/(gamma*rho(p))+
     &                       csound(nb(i))**2/(gamma*rho(nb(i)))
                 ENDIF
                 aij(i)=(gradprho+piij)*0.5*(dwmass+dwmass1)
                 muij(i)=ABS(muij(i))
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 40 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg/(dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    gradprho=2.*csound(p)*csound(nb(i))/(gamma*
     &                       SQRT(rhoprhoi))
                 ELSE
                    gradprho=csound(p)**2/(gamma*rho(p))+
     &                       csound(nb(i))**2/(gamma*rho(nb(i)))
                 ENDIF
                 aij(i)=(gradprho+piij)*dwmass
                 muij(i)=ABS(muij(i))
 40           CONTINUE

           ENDIF

           IF(isubset(p).EQ.1) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 50 i=1,inear
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 aijmass=aij(i)*mass(nb(i))*cosmofac
                 acc(p,1)=acc(p,1)-aijmass*dx
                 acc(p,2)=acc(p,2)-aijmass*dy
                 acc(p,3)=acc(p,3)-aijmass*dz
 50           CONTINUE
           ENDIF

           CALL WHENEQ(inear,jsubset,1,1,nnind,ninear)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,ninear
              j=nnind(i)
              dx=pos(p,1)-pos(nb(j),1)
              dy=pos(p,2)-pos(nb(j),2)
              dz=pos(p,3)-pos(nb(j),3)
              IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
              IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
              IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
              IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
              IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
              IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
              aijmass=aij(j)*mass(p)*cosmofac
              acc(nb(j),1)=acc(nb(j),1)+aijmass*dx
              acc(nb(j),2)=acc(nb(j),2)+aijmass*dy
              acc(nb(j),3)=acc(nb(j),3)+aijmass*dz
 60        CONTINUE

 90     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
          SUBROUTINE bisect(nblist,tempmin,tempmax,smallnum,dt)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the thermal energy implicitly by a
C     timestep dt, using a bisection method.  The local variable p
C     is a pointer to the bodies.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER i,ib,nblist,nitb,nitbmax
        LOGICAL fixpdv
        REAL dt,smallnum,tempmin,tempmax,tlast,tlow,thigh,logtlow,
     &       logthigh,tnew,logtnew,flogtl,flogth,flogtn,
     &       ncool,lcool,hcool
     
        PARAMETER(nitbmax=100)
     
C=======================================================================

        DO 20 i=1,nblist

           ib=bodlist(i)

           IF(.NOT.radiate) THEN

              dethold(ib)= -gamma1*ethold(ib)*hsmdivv(ib)/(cosmfold*
     &                     hsmooth(ib))
              dethdt(ib)=cosmfold*dethold(ib)/cosmofac
              ethermal(ib)=ethold(ib)+0.5*dt*(dethold(ib)+dethdt(ib))

              IF(ethermal(ib).LE.0.0)
     &           CALL terror(' negative thermal energy in bisect ')
C                     ------

           ELSE

              nitb=0
              fixpdv=.FALSE.
              tlow=0.05*meanmwt*mhboltz*gamma1*ethold(ib)
              thigh=20.0*meanmwt*mhboltz*gamma1*ethold(ib)
              tnew=meanmwt*mhboltz*gamma1*ethold(ib)
              IF(tlow.LT.tempmin) tlow=tempmin
              IF(thigh.GT.tempmax) thigh=tempmax
              logtlow=LOG10(tlow)
              logthigh=LOG10(thigh)
              logtnew=LOG10(tnew)

              CALL ethdiff(ib,tlow,logtlow,dt,flogtl,fixpdv,hcool)
C                  -------
              CALL ethdiff(ib,thigh,logthigh,dt,flogth,fixpdv,lcool)
C                  -------
              CALL ethdiff(ib,tnew,logtnew,dt,flogtn,fixpdv,ncool)
C                  -------
              IF(dt*abs(hcool)/ethold(ib).LT.1.e-18.AND.dt*abs(lcool)/
     &               ethold(ib).LT.1.e-18.AND.dt*abs(ncool)/ethold(ib)
     &               .LT.1.e-18.AND.dt*abs(derad(ib))/ethold(ib).LT.
     &               1.e-18) THEN
                  ethermal(ib)=ethold(ib)+0.5*dt*(dethdt(ib)+
     &                dethold(ib))
                  IF(ethermal(ib).LE.zero) THEN
                      dethold(ib)= -gamma1*ethold(ib)*hsmdivv(ib)/
     &                             (hsmooth(ib)*cosmfold)
                      dethdt(ib)=cosmfold*dethold(ib)/cosmofac
                      ethermal(ib)=ethold(ib)+0.5*dt*(dethold(ib)+
     &                       dethdt(ib))

                      IF(ethermal(ib).LE.0.0)
     &                CALL terror(' negative thermal energy in bisect ')
C                          ------
                  ENDIF
                  GO TO 15
              ENDIF

              IF(flogtl*flogth.GE.0.0) THEN

                 tnew=ethold(ib)*meanmwt*mhboltz*gamma1

                 IF(tnew.LT.tempmax/20.) THEN

                    fixpdv=.TRUE.

                    CALL ethdiff(ib,tlow,logtlow,dt,flogtl,fixpdv,
C                        -------
     &                           lcool)
                    CALL ethdiff(ib,thigh,logthigh,dt,flogth,fixpdv,
C                        -------
     &                           hcool)
                 ELSE

                    write(6,*)IB,TLOW,LOGTLOW,DT,FLOGTL,FIXPDV
                    write(6,*)IB,THIGH,LOGTHIGH,DT,FLOGTH,FIXPDV
                    write(6,*)mass(ib),rho(ib),itimestp(ib),otimestp(ib)
                    write(6,*)derad(ib),csound(ib),dethdt(ib),
     &                        dethold(ib)
                    write(6,*)ethermal(ib),mumaxdvh(ib),hsmdivv(ib)
                    write(6,*)tnow,hsmooth(ib)
                    CALL terror(' limit error 1 in bisect ')
C                        ------
                 ENDIF

              ENDIF

              tlast=tlow
     
 10           CONTINUE
     
              IF(flogtl*flogth.GE.0.0) THEN

                 IF(tnew.LT.tempmin*20.) THEN
                    ethermal(ib)=tempmin/(meanmwt*mhboltz*gamma1)
                    GO TO 20
                 ELSE

                    write(6,*)IB,TLOW,LOGTLOW,DT,FLOGTL,FIXPDV
                    write(6,*)IB,THIGH,LOGTHIGH,DT,FLOGTH,FIXPDV
                    write(6,*)mass(ib),rho(ib),itimestp(ib),otimestp(ib)
                    write(6,*)derad(ib),csound(ib),dethdt(ib),
     &                        dethold(ib)
                    write(6,*)ethermal(ib),mumaxdvh(ib),hsmdivv(ib)
                    write(6,*)tnow,hsmooth(ib)
                    CALL terror(' limit error 2 in bisect ')
C                        ------
     
                 ENDIF
              ENDIF

              logtnew=0.5*(logtlow+logthigh)
              tnew=10.**logtnew
     
              CALL ethdiff(ib,tnew,logtnew,dt,flogtn,fixpdv,ncool)
C                  -------

              IF(flogtl*flogtn.LE.0.0) THEN
                 flogth=flogtn
                 logthigh=logtnew
              ELSE
                 flogtl=flogtn
                 logtlow=logtnew
              ENDIF
     
              nitb=nitb+1
     
              IF(nitb.GT.nitbmax)
     &           CALL terror(' convergence error in bisect ')
C                     ------
     
              IF(ABS((tlast-tnew)/tlast).GT.smallnum.AND.
     &               flogtn.NE.0.)  THEN
                 tlast=tnew
                 GO TO 10
              ENDIF
     
              ethermal(ib)=tnew/(meanmwt*mhboltz*gamma1)
   
           ENDIF
 15     CONTINUE
 20     CONTINUE
     
        RETURN
        END
C***********************************************************************
C
C
          SUBROUTINE entdiff(ib,temp,logtemp,dt,flogt,fixpdv)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the difference of the entropy at the 
C     advanced time step the the sum of the entropy at the last step 
C     and the mean of the time derivatives of the entropy at both 
C     steps.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER ib
        LOGICAL fixpdv
        REAL dt,tenf1,tenf2,heatcool,tempdiff,cutoff,tendiff,tenf3,
     &       tenf3exp,temp,logtemp,flogt,compcool,aconst,tcutoff,
     &       tcutoffc,compdiff,comdiff,cutoffc
     
C=======================================================================

        IF(fixpdv) THEN
           dentold(ib)= zero
           dentdt(ib)= dentold(ib)
        ENDIF

        tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(ib)*nsmooth/8./
     &          (rho(ib)/cosmof3))**(2./3.))-((epsvect(ib)/
     &          piinv)**2.)/3.)*(rho(ib)/cosmof3)
        IF(tcutoff.LT.mintemp) tcutoff=mintemp
        tempdiff=10.0*(temp-tcutoff)/(tcutoff+tiny)
        IF(tempdiff.GT.thirty5) tempdiff=thirty5
        IF(tempdiff.LT.-thirty5) tempdiff = -thirty5
        tendiff=10.**tempdiff
        cutoff=0.5*(tendiff-1.0)/(tendiff+1.0)+0.5

        tcutoffc=tcutoff
        IF(tcutoff.LT.10000.) tcutoffc=10000.
        compdiff=10.0*(temp-tcutoffc)/(tcutoffc+tiny)
        IF(compdiff.GT.thirty5) compdiff=thirty5
        IF(compdiff.LT.-thirty5) compdiff=-thirty5
        comdiff=10.**compdiff
        cutoffc=0.5*(comdiff-1.0)/(comdiff+1.0)+0.5

        tenf1=10.**(-4.43-0.273*(4.-logtemp)**2)
        tenf2=10.**(-0.10-1.880*(5.23-logtemp)**4)
        tenf3exp=6.2-logtemp
        IF(logtemp-6.2.GT.zero) tenf3exp=zero
        tenf3=10.**(-2.*tenf3exp**4-1.7)
        compcool=comptmh*temp*((1.+redshift)**4)

C##### LH 050990

        heatcool=aheatmh*fhydrogn+bheatmh2*rho(ib)*fhydrog2/
     &           cosmof3-ccoolmh2*rho(ib)*(tenf1+tenf2+tenf3)*
     &           cutoff*fhydrog2/cosmof3-compcool*cutoffc

C##### LH 050990

        aconst=slowcool*temp/(meanmwt*mhboltz*gamma1*dt)+
     &         rho(ib)**gamma1*dentdt(ib)/gamma1
        IF(aconst.LE.zero) aconst=tiny
        aconst=heatcool/aconst
        heatcool=heatcool/SQRT(1.+aconst*aconst)
        flogt=temp/(meanmwt*mhboltz*rho(ib)**gamma1)-entold(ib)-
     &        0.5*dt*(dentdt(ib)+gamma1*heatcool/rho(ib)**gamma1+
     &        dentold(ib)+gamma1*derad(ib)/rho(ib)**gamma1)

        RETURN
        END
C***********************************************************************
C
C
         SUBROUTINE bisectas(nblist,tempmin,tempmax,smallnum,dt)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the entropic function a(s) implicitly by 
C     a timestep dt, using a bisection method.  The local variable p
C     is a pointer to the bodies.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER i,ib,nblist,nitb,nitbmax
        LOGICAL fixpdv
        REAL dt,smallnum,tempmin,tempmax,tlast,tlow,thigh,logtlow,
     &       logthigh,tnew,logtnew,flogtl,flogth,flogtn
     
        PARAMETER(nitbmax=100)
     
C=======================================================================

        DO 20 i=1,nblist

           ib=bodlist(i)

           IF(.NOT.radiate) THEN

              dentold(ib)= zero
              dentdt(ib)=dentold(ib)
              entropy(ib)=entold(ib)+0.5*dt*(dentold(ib)+dentdt(ib))

              IF(entropy(ib).LE.0.0)
     &           CALL terror(' negative entropy in bisectas ')
C                     ------

           ELSE

              nitb=0
              fixpdv=.FALSE.
              tlow=0.1*entropy(ib)*meanmwt*mhboltz*rho(ib)**gamma1
              thigh=20.0*entropy(ib)*meanmwt*mhboltz*rho(ib)**gamma1
              IF(tlow.LT.tempmin) tlow=tempmin
              IF(thigh.GT.tempmax) thigh=tempmax
              logtlow=LOG10(tlow)
              logthigh=LOG10(thigh)

              CALL entdiff(ib,tlow,logtlow,dt,flogtl,fixpdv)
C                  -------
              CALL entdiff(ib,thigh,logthigh,dt,flogth,fixpdv)
C                  -------

              IF(flogtl*flogth.GE.0.0) THEN

                 tnew=entropy(ib)*meanmwt*mhboltz*rho(ib)**gamma1

                 IF(tnew.LT.tempmin*10.) THEN

                    fixpdv=.TRUE.

                    CALL entdiff(ib,tlow,logtlow,dt,flogtl,fixpdv)
C                        -------
                    CALL entdiff(ib,thigh,logthigh,dt,flogth,fixpdv)
C                        -------
                 ELSE

                    CALL terror(' limit error 1 in bisectas ')
C                        ------
                 ENDIF

              ENDIF

              tlast=tlow
     
 10           CONTINUE
     
              IF(flogtl*flogth.GE.0.0) THEN

                 IF(tnew.LT.tempmin*10.) THEN
                    entropy(ib)=tempmin/(meanmwt*mhboltz*
     &                          rho(ib)**gamma1)
                    GO TO 20
                 ELSE

                    CALL terror(' limit error 2 in bisectas ')
C                        ------

                 ENDIF
              ENDIF
     
              logtnew=0.5*(logtlow+logthigh)
              tnew=10.**logtnew
     
              CALL entdiff(ib,tnew,logtnew,dt,flogtn,fixpdv)
C                  -------

              IF(flogtl*flogtn.LE.0.0) THEN
                 flogth=flogtn
                 logthigh=logtnew
              ELSE
                 flogtl=flogtn
                 logtlow=logtnew
              ENDIF
     
              nitb=nitb+1
     
              IF(nitb.GT.nitbmax)
     &           CALL terror(' convergence error in bisectas ')
C                     ------
     
              IF(ABS((tlast-tnew)/tlast).GT.smallnum)  THEN
                 tlast=tnew
                 GO TO 10
              ENDIF
     
              entropy(ib)=tnew/(meanmwt*mhboltz*rho(ib)**gamma1)
   
           ENDIF
  
 20     CONTINUE
     
        RETURN
        END
C***********************************************************************
C
C
          SUBROUTINE ethdiff(ib,temp,logtemp,dt,flogt,fixpdv,heatcool)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the difference of the thermal energy at
C     the advanced time step the the sum of the thermal energy at 
C     the last step and the mean of the time derivatives of the
C     thermal at both steps.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER ib
        LOGICAL fixpdv
        REAL dt,tenf1,tenf2,heatcool,tempdiff,cutoff,tenf3,
     &       tenf3exp,temp,logtemp,flogt,compcool,aconst,tcutoff,
     &       tcutoffc,compdiff,cutoffc
     
C=======================================================================

        IF(fixpdv) THEN
           dethold(ib)= -gamma1*ethold(ib)*hsmdivv(ib)/(cosmofac*
     &                  hsmooth(ib))
           dethdt(ib)= dethold(ib)
        ENDIF

        tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(ib)*nsmooth/8./
     &          (rho(ib)/cosmof3))**(2./3.))-((epsvect(ib)/
     &          piinv)**2.)/3.)*(rho(ib)/cosmof3)
        IF(tcutoff.LT.mintemp) tcutoff=mintemp
        tempdiff= 30.*(1.-temp/tcutoff)
        IF(tempdiff.LT.0.)tempdiff=0.
        cutoff=exp(-tempdiff *tempdiff)

        tcutoffc=tcutoff
        IF(tcutoff.LT.10000.) tcutoffc=10000.
        compdiff= 30.*(1.-temp/tcutoffc)
        IF(compdiff.LE.0.)compdiff=0.
        cutoffc=exp(-compdiff *compdiff)

        tenf1=10.**(-4.43-0.273*(4.-logtemp)**2)
        tenf2=10.**(-0.10-1.880*(5.23-logtemp)**4)
        tenf3exp=6.2-logtemp
        IF(logtemp-6.2.GT.zero) tenf3exp=zero
        tenf3=10.**(-2.*tenf3exp**4-1.7)
        compcool=comptmh*temp*((1.+redshift)**4)
        heatcool=aheatmh*fhydrogn+bheatmh2*rho(ib)*fhydrog2/
     &           cosmof3-ccoolmh2*rho(ib)*(tenf1+tenf2+tenf3)*
     &           cutoff*fhydrog2/cosmof3-compcool*cutoffc
        aconst=slowcool*temp/(meanmwt*mhboltz*gamma1*dt)+dethdt(ib)
        IF(aconst.LE.zero) aconst=tiny
        aconst=heatcool/aconst
        heatcool=heatcool/SQRT(1.+aconst*aconst)
        flogt=temp/(meanmwt*mhboltz*gamma1)-ethold(ib)-0.5*dt*
     &        (dethdt(ib)+heatcool+dethold(ib)+derad(ib))

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE gravity(option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute gravitational potential and acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option

C=======================================================================

        IF(selfgrav) THEN

           IF(variabls.AND.nsph.GT.0) THEN

              CALL maketeps
C                  --------
              CALL stepeps
C                  -------
           ENDIF

           IF((.NOT.dgravsum).OR.variabls) CALL maketree
C                                               --------

           IF(variabls.AND.nsph.EQ.0) THEN

              CALL maketeps
C                  --------
              CALL stepeps
C                  -------
           ENDIF

           IF(.NOT.dgravsum) CALL celleps
C                                 -------
           CALL accgrav(option)
C               -------

        ENDIF


        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE maketeps
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the tree structure for computing
C     variable gravitational softening lengths for collisionless
C     particles.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

        LOGICAL firstc

        SAVE firstc

        DATA firstc/.TRUE./

C=======================================================================

        IF(nsph.GT.0.OR.firstc) THEN

           firstc=.FALSE.

           CALL setbox('coll')
C               ------
           CALL loadtree('coll')
C               --------
           CALL cellnumb
C               --------
        ENDIF

C   Compute coordinates of edges of cells.
C   --------------------------------------
        CALL celledge
C            --------

C   Determine which cells contain ingroup particles.
C   ------------------------------------------------
        CALL groupcel
C            --------

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stepeps
C
C
C***********************************************************************
C
C
C     Subroutine to update softening lengths for all collisionless
C     particles.  The list of neighbors is returned from the 
C     subroutine findnear in the vector nearlist, which is 
C     equivalenced to the common array bodlist.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,pbody,npnear,inear,nnearbod,nearpb(nbodsmax),
     &          np,jsubset(nbodsmax),nearlist(nbodsmax),
     &          keepnear(nbodsmax),nnearfix
        REAL testnn,foureps2,dx,dy,dz

        EQUIVALENCE (nearlist(1),bodlist(1)),(jsubset(1),subindex(1)),
     &              (nearpb(1),parent(1)),(keepnear(1),asubp(1))

C=======================================================================

C   Extrapolate softening lengths.
C   ------------------------------
        DO 10 p=nsph+1,nbodies
           IF(nnear(p).EQ.0) THEN
              nnearfix=1
           ELSE
              nnearfix=0
           ENDIF
           epsvect(p)=epsvect(p)*0.5*((REAL(nsvolume)/REAL(nnear(p)+
     &                nnearfix))**(one/3.)+one)
 10     CONTINUE

C   Initialize neighbor diagnostics.
C   --------------------------------
        nstot=0
        nsmin=nbodies
        nsmax=0

C   Find nearest neighbors of grouped cells.
C   ----------------------------------------

        DO 100 p=1,ngroups

C   Find nearest neighbors.
C   -----------------------
           CALL findnear(p,npnear)
C               --------

           DO 70 np=pgroupb(p),pgroupb(p+1)-1

              pbody=groupbod(np)

              foureps2=four*epsvect(pbody)*epsvect(pbody)

              DO 20 i=1,npnear
                 dx=pos(pbody,1)-pos(nearlist(i),1)
                 dy=pos(pbody,2)-pos(nearlist(i),2)
                 dz=pos(pbody,3)-pos(nearlist(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 tempvect(i)=dx**2+dy**2+dz**2
                 IF(nearlist(i).EQ.pbody) tempvect(i)=foureps2+one
 20           CONTINUE

C   Filter out neighbors further than two smoothing lengths.
C   --------------------------------------------------------

              CALL WHENFLT(npnear,tempvect,1,foureps2,isubset,inear)

              testnn=ABS(REAL(inear-nsvolume)/REAL(nsvolume))

              IF(testnn.LE.nsvtol) THEN

                 nstot=nstot+inear
                 IF(inear.LT.nsmin) nsmin=inear
                 IF(inear.GT.nsmax) nsmax=inear
                 nnear(pbody)=inear

              ELSE

                 IF(inear.GT.nsvolume) THEN

                    DO 50 i=1,inear
                       nearpb(i)=nearlist(isubset(i))
 50                 CONTINUE
 
                    CALL reducee(pbody,inear,nnearbod)
C                        -------

                 ELSE

                    CALL enlargee(pbody,nnearbod)
C                        --------
                 ENDIF

              ENDIF

 70        CONTINUE

 100    CONTINUE

        IF(nbodies.NE.nsph) THEN
           nsavg=nstot/(nbodies-nsph)
        ELSE
           nsavg=0
        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                      SUBROUTINE findnear(p,npnear)
C
C
C***********************************************************************
C
C
C     Subroutine to search for nearest neighbors of the grouped
C     cell p.  Vectorization is achieved by processing all cells
C     at a given level in the tree simultaneously.  The list of
C     neighbors is returned through the vector nearlist, which is
C     equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nsubdiv,nearlist(nbodsmax),npnear,nnodes,nkeep,
     &          nodelist(ncells),keepnear(nbodsmax),keepstak(nbodsmax)
        LOGICAL testcrit
        REAL pbottom(ndim),ptop(ndim),dx,dy,dz,xnode,ynode,znode,
     &       groupx,groupy,groupz

        EQUIVALENCE (nearlist(1),bodlist(1)),(nodelist(1),celllist(1)),
     &              (keepnear(1),parent(1)),(keepstak(1),asubp(1))

C=======================================================================

        npnear=0

        CALL srchbox(p,npnear,pbottom,ptop)
C            -------

        nnodes=intone
        nodelist(1)=root

 20     CONTINUE

        IF(nnodes.GT.0) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nnodes
              xnode=bottom(nodelist(i),1)+0.5*cellsize(nodelist(i))
              ynode=bottom(nodelist(i),2)+0.5*cellsize(nodelist(i))
              znode=bottom(nodelist(i),3)+0.5*cellsize(nodelist(i))
              groupx=0.5*(ptop(1)+pbottom(1))
              groupy=0.5*(ptop(2)+pbottom(2))
              groupz=0.5*(ptop(3)+pbottom(3))
              dx=0.0
              dy=0.0
              dz=0.0
              IF(xnode-groupx.GE.hboxsize.AND.bwrap) dx=-pboxsize
              IF(xnode-groupx.LT.-hboxsize.AND.bwrap) dx=pboxsize
              IF(ynode-groupy.GE.hboxsize.AND.bwrap) dy=-pboxsize
              IF(ynode-groupy.LT.-hboxsize.AND.bwrap) dy=pboxsize
              IF(znode-groupz.GE.hboxsize.AND.bwrap) dz=-pboxsize
              IF(znode-groupz.LT.-hboxsize.AND.bwrap) dz=pboxsize
              testcrit=     pbottom(1).LE.(bottom(nodelist(i),1)+
     &                                   cellsize(nodelist(i))+dx)
     &                 .AND.pbottom(2).LE.(bottom(nodelist(i),2)+
     &                                   cellsize(nodelist(i))+dy)
     &                 .AND.pbottom(3).LE.(bottom(nodelist(i),3)+
     &                                   cellsize(nodelist(i))+dz)
     &                 .AND.ptop(1).GE.(bottom(nodelist(i),1)+dx)
     &                 .AND.ptop(2).GE.(bottom(nodelist(i),2)+dy)
     &                 .AND.ptop(3).GE.(bottom(nodelist(i),3)+dz)
     &                 .AND.nodelist(i).NE.groups(p)
              IF(testcrit.AND.nodelist(i).LE.nbodsmax) THEN
                 keepnear(i)=2
              ELSE
                 keepnear(i)= -2
              ENDIF
              IF(testcrit.AND.nodelist(i).GT.nbodsmax) THEN
                 keepstak(i)=2
              ELSE
                 keepstak(i)= -2
              ENDIF
 30        CONTINUE

           CALL WHENIGT(nnodes,keepnear,1,0,isubset,nkeep)

           IF(npnear+nkeep.GT.nbodsmax)
     &        CALL terror(' array overflow in findnear ')
C                  ------

           DO 40 i=1,nkeep
              nearlist(npnear+i)=nodelist(isubset(i))
 40        CONTINUE

           npnear=npnear+nkeep

           CALL WHENIGT(nnodes,keepstak,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.npart)THEN
	      WRITE(uterm,*)8*nsubdiv
              CALL terror(' asubp overflow in findnear ')
C                  ------
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 50        CONTINUE

           CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

           IF(nnodes.GT.nbodsmax.OR.nnodes.GT.ncells)THEN
	      WRITE(uterm,*)8*nsubdiv,nnodes
              CALL terror(' nodelist overflow in findnear ')
C                  ------
           ENDIF

           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
               SUBROUTINE reducee(pbody,inear,nnearbod)
C
C
C***********************************************************************
C
C
C     Subroutine to decrease softening length of body pbody so that
C     the softening volume contains roughly nsvolume neighbors.  The
C     list of nearest neighbors found in stepeps is passed through
C     the vector nearpb, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,pbody,inear,nhtest,nearpb(nbodsmax),nnearbod,
     &          jsubset(nbodsmax),keepnear(nbodsmax),nhmin,nhmax
        REAL hmin,hmax,htest,fourhsm2,dx,dy,dz

        EQUIVALENCE (nearpb(1),parent(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1))

C=======================================================================

        hmax=epsvect(pbody)

        nhmax=inear

        DO 5 i=1,inear
           isubset(i)=i
           dx=pos(pbody,1)-pos(nearpb(i),1)
           dy=pos(pbody,2)-pos(nearpb(i),2)
           dz=pos(pbody,3)-pos(nearpb(i),3)
           IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
           IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
           IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
           IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
           IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
           IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
           tempvect(i)=dx**2+dy**2+dz**2
 5      CONTINUE

        hmin=0.5*hmax

 7      CONTINUE

        htest=hmin

        fourhsm2=four*htest*htest

        CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)

        IF(nhtest.GE.nsvolume) THEN
           hmin=0.5*hmin
           GO TO 7
        ENDIF

        nhmin=nhtest

 10     CONTINUE

        IF(ABS(REAL(nhtest-nsvolume)).GT.nsvtol*REAL(nsvolume)) THEN

           htest=0.5*(hmin+hmax)
           fourhsm2=four*htest*htest

           CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)

           IF(nhmin.EQ.nhmax)
     &        CALL terror(' iteration error in reducee ')
C                  ------

           IF(nhtest.GT.nsvolume) THEN
              hmax=htest
              nhmax=nhtest
           ELSE
              hmin=htest
              nhmin=nhtest
           ENDIF

           GO TO 10

        ENDIF

        nstot=nstot+nhtest
        IF(nhtest.LT.nsmin) nsmin=nhtest
        IF(nhtest.GT.nsmax) nsmax=nhtest
        nnear(pbody)=nhtest
        epsvect(pbody)=htest

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE enlargee(pbody,nnearbod)
C
C
C***********************************************************************
C
C
C     Subroutine to increase softening length of body pbody so that
C     the softening volume contains roughly nsvolume neighbors.  The
C     neighbor list is returned from the subroutine nearpeps in the
C     vector nearlist, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,pbody,nhtest,nearlist(nbodsmax),nnearbod,
     &          jsubset(nbodsmax),npnear,keepnear(nbodsmax),nhold
        REAL hmin,hmax,htest,fourhsm2,testnn,hsearch,dx,dy,dz

        EQUIVALENCE (nearlist(1),parent(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1))

C=======================================================================

        hsearch=2.0
        npnear=0

 5      CONTINUE

        IF(npnear.LT.nsvolume) THEN

           hsearch=1.5*hsearch

           IF(hsearch.GT.1.e5)
     &        CALL terror(' search error in enlargee ')
C                  ------

           CALL nearpeps(hsearch,pbody,npnear)
C               --------

           GO TO 5

        ENDIF
        
        nhtest=npnear

        DO 10 i=1,npnear
           dx=pos(pbody,1)-pos(nearlist(i),1)
           dy=pos(pbody,2)-pos(nearlist(i),2)
           dz=pos(pbody,3)-pos(nearlist(i),3)
           IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
           IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
           IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
           IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
           IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
           IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
           tempvect(i)=dx**2+dy**2+dz**2
           isubset(i)=i
 10     CONTINUE

        testnn=ABS(REAL(npnear-nsvolume))/REAL(nsvolume)

        IF(npnear.LE.nsvolume.OR.testnn.LE.nsvtol) THEN

           htest=hsearch*epsvect(pbody)

        ELSE

           hmin=0.
           hmax=hsearch*epsvect(pbody)

           nhold=-1

 20        CONTINUE

           IF(ABS(REAL(nhtest-nsvolume)).GT.nsvtol*REAL(nsvolume)) THEN

              htest=0.5*(hmin+hmax)
              fourhsm2=four*htest*htest

              CALL WHENFLT(npnear,tempvect,1,fourhsm2,isubset,nhtest)

C              IF(nhtest.EQ.nhold)
C     &           CALL terror(' iteration error in reducee ')
C                     ------

              nhold=nhtest

              IF(nhtest.GT.nsvolume) THEN
                 hmax=htest
              ELSE
                 hmin=htest
              ENDIF

              GO TO 20

           ENDIF

        ENDIF

        nstot=nstot+nhtest
        IF(nhtest.LT.nsmin) nsmin=nhtest
        IF(nhtest.GT.nsmax) nsmax=nhtest
        nnear(pbody)=nhtest
        epsvect(pbody)=htest

        RETURN
        END
C***********************************************************************
C
C
                SUBROUTINE nearpeps(hsearch,pbody,npnear)
C
C
C***********************************************************************
C
C
C     Subroutine to search for nearest neighbors of body pbody within
C     hsearch softening lengths.  Vectorization is achieved by 
C     processing all cells at a given level in the tree simultaneously.
C     The neighbor list is passed back to the calling routine in the
C     vector nearlist, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER pbody,i,nsubdiv,nearlist(nbodsmax),npnear,nnodes,nkeep,
     &          k,ibody,jsubset(nbodsmax),njsubset,nodelist(ncells),
     &          keepnear(nbodsmax),keepstak(nbodsmax)
        LOGICAL testcrit
        REAL pbottom(ndim),ptop(ndim),hsearch,sradius,dx,dy,dz,
     &       xnode,ynode,znode

        EQUIVALENCE (nearlist(1),parent(1)),(jsubset(1),subindex(1)),
     &              (nodelist(1),celllist(1)),(keepnear(1),subindex(1)),
     &              (keepstak(1),asubp(1))

C=======================================================================

        npnear=0
        sradius=hsearch*epsvect(pbody)

        DO 10 k=1,3
           ptop(k)=pos(pbody,k)+sradius
           pbottom(k)=pos(pbody,k)-sradius
 10     CONTINUE

        nnodes=intone
        nodelist(1)=root

 20     CONTINUE

        IF(nnodes.GT.0) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nnodes
              xnode=bottom(nodelist(i),1)+0.5*cellsize(nodelist(i))
              ynode=bottom(nodelist(i),2)+0.5*cellsize(nodelist(i))
              znode=bottom(nodelist(i),3)+0.5*cellsize(nodelist(i))
              dx=0.0
              dy=0.0
              dz=0.0
              IF(xnode-pos(pbody,1).GE.hboxsize.AND.bwrap) dx=-pboxsize
              IF(xnode-pos(pbody,1).LT.-hboxsize.AND.bwrap) dx=pboxsize
              IF(ynode-pos(pbody,2).GE.hboxsize.AND.bwrap) dy=-pboxsize
              IF(ynode-pos(pbody,2).LT.-hboxsize.AND.bwrap) dy=pboxsize
              IF(znode-pos(pbody,3).GE.hboxsize.AND.bwrap) dz=-pboxsize
              IF(znode-pos(pbody,3).LT.-hboxsize.AND.bwrap) dz=pboxsize
              testcrit=     pbottom(1).LE.(bottom(nodelist(i),1)+
     &                                   cellsize(nodelist(i))+dx)
     &                 .AND.pbottom(2).LE.(bottom(nodelist(i),2)+
     &                                   cellsize(nodelist(i))+dy)
     &                 .AND.pbottom(3).LE.(bottom(nodelist(i),3)+
     &                                   cellsize(nodelist(i))+dz)
     &                 .AND.ptop(1).GE.(bottom(nodelist(i),1)+dx)
     &                 .AND.ptop(2).GE.(bottom(nodelist(i),2)+dy)
     &                 .AND.ptop(3).GE.(bottom(nodelist(i),3)+dz)
     &                 .AND.nodelist(i).NE.pbody
              IF(testcrit.AND.nodelist(i).LE.nbodsmax) THEN
                 keepnear(i)=2
              ELSE
                 keepnear(i)= -2
              ENDIF
              IF(testcrit.AND.nodelist(i).GT.nbodsmax) THEN
                 keepstak(i)=2
              ELSE
                 keepstak(i)= -2
              ENDIF
 30        CONTINUE

           CALL WHENIGT(nnodes,keepnear,1,0,isubset,nkeep)

           DO 35 i=1,nkeep
              ibody=nodelist(isubset(i))
              dx=pos(pbody,1)-pos(ibody,1)
              dy=pos(pbody,2)-pos(ibody,2)
              dz=pos(pbody,3)-pos(ibody,3)
              IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
              IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
              IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
              IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
              IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
              IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
              tempvect(i)=dx**2+dy**2+dz**2
 35        CONTINUE

           CALL WHENFLT(nkeep,tempvect,1,sradius**2,jsubset,njsubset)

           nkeep=njsubset

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nkeep
              nearlist(npnear+i)=nodelist(isubset(jsubset(i)))
 40        CONTINUE

           npnear=npnear+nkeep

           CALL WHENIGT(nnodes,keepstak,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.nbodsmax.OR.8*nsubdiv.GT.ncells)
     &        CALL terror(' asubp overflow in nearpeps ')
C                  ------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 50        CONTINUE

           CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE maketree
C
C
C***********************************************************************
C
C
C     Main routine to control initialization of the tree structure 
C     for computing the gravitational interaction.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================

C   Set box properties.
C   -------------------
        CALL setbox('all ')
C            ------
 
C   Load bodies into the tree.
C   --------------------------
        CALL loadtree('all ')
C            --------

C   Compute properties of cells.
C   ----------------------------
        CALL hackcell
C            --------
 
        incellsg=incells
        ttree=tpos

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE hackcell
C
C
C***********************************************************************
C
C
C     Subroutine to compute masses, center of mass coordinates,
C     and optional quadrupole moments of cells, processing cells
C     in order of increasing size.  The permutation vector is
C     stored in the common variable celllist.  Vectorization is
C     achieved by simultaneously processing all cells at the
C     same level in the hierarchy.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,k,l,m,n,nsubb,nsubc,nnodes,mupper,
     &          lcf1

C=======================================================================
        
C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize properties of cells.
C   -------------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           mass(p)=0.
           npercell(p)=0
           pos(p,1)=0.
           pos(p,2)=0.
           pos(p,3)=0.
 10     CONTINUE

        IF(usequad) THEN
           DO 30 k=1,2*ndim-1
              DO 20 p=nbodsmax+1,nbodsmax+incells
                 quad(p,k)=0.
 20           CONTINUE
 30        CONTINUE
        ENDIF

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 40     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 50 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 60
              ENDIF
 50        CONTINUE                    

 60        CONTINUE

           lcf1=lcell-fcell+1

           IF(lcf1.GT.ncells.OR.lcf1.GT.nbodsmax)
     &        CALL terror(' lcf1 overflow in hackcell ')
C                  ------

C   Compute properties of the selected cells, looping over subcells.
C   ----------------------------------------------------------------

           DO 110 j=1,nsubcell

              DO 70 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

              IF(nnodes.GT.ncells.OR.nnodes.GT.nbodsmax)
     &           CALL terror(' array overflow in hackcell ')
C                     ------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
                 mass(parent(i))=mass(parent(i))+mass(asubp(i))
                 pos(parent(i),1)=pos(parent(i),1)+mass(asubp(i))*
     &                            pos(asubp(i),1)
                 pos(parent(i),2)=pos(parent(i),2)+mass(asubp(i))*
     &                            pos(asubp(i),2)
                 pos(parent(i),3)=pos(parent(i),3)+mass(asubp(i))*
     &                            pos(asubp(i),3)
 80           CONTINUE

              CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 90 i=1,nsubc
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+
     &              npercell(asubp(isubset(i)))
 90           CONTINUE

              CALL WHENILE(nnodes,asubp,1,nbodsmax,isubset,nsubb)

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 100 i=1,nsubb
                 templist(i)=parent(isubset(i))
                 npercell(templist(i))=npercell(templist(i))+1
 100          CONTINUE

 110       CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 120 i=fcell,lcell
              pos(celllist(i),1)=pos(celllist(i),1)/mass(celllist(i))
              pos(celllist(i),2)=pos(celllist(i),2)/mass(celllist(i))
              pos(celllist(i),3)=pos(celllist(i),3)/mass(celllist(i))
 120       CONTINUE

C   Compute optional quadrupole moments.
C   ------------------------------------

           IF(usequad) THEN

              DO 210 j=1,nsubcell

                 DO 130 i=fcell,lcell
                    asubp(i-fcell+1)=subp(celllist(i),j)
 130             CONTINUE

                 CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 140 i=1,nnodes
                    parent(i)=celllist(isubset(i)+fcell-1)
                    asubp(i)=subp(parent(i),j)
 140             CONTINUE

                 CALL WHENIGT(nnodes,asubp,1,nbodsmax,isubset,nsubc)

                 IF(ndim.GT.2) THEN
                    mupper=2
                 ELSE
                    mupper=ndim
                 ENDIF

                 DO 200 m=1,mupper
                    DO 190 n=m,ndim

                       l=(m-1)*(ndim-1)+n

CVD$ NODEPCHK
CDIR$ IVDEP
                       DO 150 i=1,nnodes
                          quad(parent(i),l)=quad(parent(i),l)+
     &                       mass(asubp(i))*(3.*(pos(asubp(i),m)-
     &                       pos(parent(i),m))*(pos(asubp(i),n)-
     &                       pos(parent(i),n)))
 150                   CONTINUE

                       IF(m.EQ.n) THEN
                          DO 170 k=1,ndim
CVD$ NODEPCHK
CDIR$ IVDEP
                             DO 160 i=1,nnodes
                                quad(parent(i),l)=quad(parent(i),l)-
     &                             mass(asubp(i))*(pos(asubp(i),k)-
     &                             pos(parent(i),k))**2
 160                         CONTINUE
 170                      CONTINUE
                       ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
                       DO 180 i=1,nsubc
                          templist(i)=parent(isubset(i))
                          quad(templist(i),l)=quad(templist(i),l)+
     &                       quad(asubp(isubset(i)),l)
 180                   CONTINUE

 190                CONTINUE
 200             CONTINUE

 210          CONTINUE

           ENDIF

           fcell=lcell+1
             
           GO TO 40

        ENDIF

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE neighbor(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to establish nearest neighbor lists for all SPH
C     particles.  The list of neighbors is returned from the
C     subroutine findnear in the vector nearlist, which is 
C     equivalenced to the common array bodlist.  If argument pc is
C     'predict' then this subroutine performs a simple neighbor
C     search for a given set of smoothing lengths.  If pc is 'correct'
C     this subroutine will additionally adjust smoothing lengths so
C     that the number of neighbors is approximately nsmooth.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,pbody,npnear,np,inear,insave,jsubset(nbodsmax),
     &          nnearbod,nearlist(nbodsmax),keepnear(nbodsmax),
     &          nearpb(nbodsmax)
        LOGICAL keepcrit,savecrit
        REAL hsminv,testnn,dx,dy,dz

        EQUIVALENCE (nearlist(1),bodlist(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1)),(nearpb(1),parent(1))

C=======================================================================

        IF(pc.EQ.'correct') THEN

           DO 10 p=1,nsph
              templist(p)=intone
 10        CONTINUE

        ENDIF

C   Initialize nearest neighbor diagnostics.
C   ----------------------------------------
        nnearbod=0
        nntot=0
        nnmin=nsphmax
        nnmax=0

C   Find nearest neighbors of grouped cells.
C   ----------------------------------------

        DO 100 p=1,ngroups

C   Find nearest neighbors.
C   -----------------------
           CALL findnear(p,npnear)
C               --------

           DO 70 np=pgroupb(p),pgroupb(p+1)-1

              pbody=groupbod(np)
              hsminv=1./hsmooth(pbody)

              DO 20 i=1,npnear
                 dx=pos(pbody,1)-pos(nearlist(i),1)
                 dy=pos(pbody,2)-pos(nearlist(i),2)
                 dz=pos(pbody,3)-pos(nearlist(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 tempvect(i)=(dx**2+dy**2+dz**2)*hsminv*hsminv
                 IF(nearlist(i).EQ.pbody) tempvect(i)=5.0
 20           CONTINUE

C   Filter out neighbors further than two smoothing lengths.
C   --------------------------------------------------------
              CALL WHENFLT(npnear,tempvect,1,four,isubset,inear)

              testnn=ABS(REAL(inear-nsmooth)/REAL(nsmooth))
              savecrit=(.NOT.variablh).OR.(testnn.LE.nsmtol).OR.
     &                 (pc.EQ.'predict')

              IF(savecrit) THEN

                 IF(pc.EQ.'correct') templist(pbody)=-1

                 nntot=nntot+inear
                 IF(inear.LT.nnmin) nnmin=inear
                 IF(inear.GT.nnmax) nnmax=inear
                 nnear(pbody)=inear

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 30 i=1,inear
                    jsubset(i)=nearlist(isubset(i))
                    keepcrit=hsmooth(pbody).GT.hsmooth(jsubset(i)).OR.
     &                       ((hsmooth(pbody).EQ.hsmooth(jsubset(i)))
     &                       .AND.pbody.LT.jsubset(i))
                    IF(keepcrit) THEN
                       keepnear(i)=2
                    ELSE
                       keepnear(i)= -2
                    ENDIF
 30              CONTINUE

                 CALL WHENIGT(inear,keepnear,1,0,jsubset,insave)

                 pnear(pbody)=nnearbod+1
                 nnearlis(pbody)=insave
                 nnearbod=nnearbod+insave

                 IF(nnearbod.GT.maxnearb)
     &              CALL terror(' array overflow in neighbor ')
C                        ------

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 40 i=1,insave
                    nearbods(nnearbod-insave+i)=
     &                       nearlist(isubset(jsubset(i)))
 40              CONTINUE

              ELSE

                 IF(inear.GT.nsmooth) THEN

                    DO 50 i=1,inear
                       nearpb(i)=nearlist(isubset(i))
 50                 CONTINUE

                    CALL reduceh(pbody,inear,nnearbod)
C                        -------

                 ELSE

                    CALL enlargeh(pbody,nnearbod)
C                        --------
                 ENDIF

              ENDIF

 70        CONTINUE

 100    CONTINUE

        nnavg=nntot/nsph

        IF(variablg.AND.(pc.EQ.'correct')) THEN

           DO 200 i=1,nsph
              epsvect(i)=hsmooth(i)
 200       CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
               SUBROUTINE reduceh(pbody,inear,nnearbod)
C
C
C***********************************************************************
C
C
C     Subroutine to decrease smoothing length of body pbody so that
C     the smoothing volume contains roughly nsmooth neighbors.  The
C     list of nearest neighbors found in stepnear is passed through
C     the vector nearpb, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,pbody,inear,nhtest,ncheck,nearpb(nbodsmax),nnearbod,
     &          insave,jsubset(nbodsmax),keepnear(nbodsmax),nhmin,nhmax
        LOGICAL keepcrit,keepcold
        REAL hmin,hmax,htest,hsmc,dt,hsmold,fourhsm2,dx,dy,dz

        EQUIVALENCE (nearpb(1),parent(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1))

C=======================================================================

        hmax=hsmooth(pbody)

        nhmax=inear

        DO 5 i=1,inear
           isubset(i)=i
           dx=pos(pbody,1)-pos(nearpb(i),1)
           dy=pos(pbody,2)-pos(nearpb(i),2)
           dz=pos(pbody,3)-pos(nearpb(i),3)
           IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
           IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
           IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
           IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
           IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
           IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
           tempvect(i)=dx**2+dy**2+dz**2
 5      CONTINUE

        hmin=0.5*hmax

 7      CONTINUE

        htest=hmin

        fourhsm2=four*htest*htest

        CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)

        IF(nhtest.GE.nsmooth) THEN
           hmin=0.5*hmin
           GO TO 7
        ENDIF

        nhmin=nhtest

 10     CONTINUE

        IF(ABS(REAL(nhtest-nsmooth)).GT.nsmtol*REAL(nsmooth)) THEN

           htest=0.5*(hmin+hmax)
           fourhsm2=four*htest*htest

           CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)

           IF(nhmin.EQ.nhmax)
     &        CALL terror(' iteration error in reduceh ')
C                  ------

           IF(nhtest.GT.nsmooth) THEN
              hmax=htest
              nhmax=nhtest
           ELSE
              hmin=htest
              nhmin=nhtest
           ENDIF

           GO TO 10

        ENDIF

        dt=tnow-teth

        hsmc=dt*(csound(pbody)+1.2*(alpha*csound(pbody)+beta*
     &       mumaxdvh(pbody)))/courant

        IF(hsmc.GT.htest) THEN
           htest=hsmc
           IF(hsmooth(pbody).LT.htest) htest=hsmooth(pbody)
           fourhsm2=four*htest*htest
           CALL WHENFLT(inear,tempvect,1,fourhsm2,isubset,nhtest)
        ENDIF

        nntot=nntot+nhtest
        IF(nhtest.LT.nnmin) nnmin=nhtest
        IF(nhtest.GT.nnmax) nnmax=nhtest
        nnear(pbody)=nhtest
        hsmold=hsmooth(pbody)
        hsmooth(pbody)=htest
        templist(pbody)=-1

C-----------------------------------------------------------------------
C   Select bodies to add to neighbor list.  Apply standard selection
C   criteria to bodies within 2*hsmooth.  Also, add bodies within the
C   range 2*hsmooth --> 2*hsmold if they are neighbors, have been
C   processed before, and failed previous selection criteria.
C-----------------------------------------------------------------------
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nhtest
           jsubset(i)=nearpb(isubset(i))
           keepcold=hsmooth(jsubset(i)).GT.hsmold.OR.
     &                 ((hsmooth(jsubset(i)).EQ.hsmold).AND.
     &                 jsubset(i).LT.pbody)
           keepcrit=((hsmooth(pbody).GT.hsmooth(jsubset(i)).OR.
     &                 ((hsmooth(pbody).EQ.hsmooth(jsubset(i))).AND.
     &                 pbody.LT.jsubset(i))).AND.templist(jsubset(i))
     &                 .GT.0).OR.((templist(jsubset(i)).LT.0).AND.
     &                 (.NOT.keepcold))
           IF(keepcrit) THEN
              keepnear(i)=2
           ELSE
              keepnear(i)= -2
           ENDIF
 30     CONTINUE

        CALL WHENIGT(nhtest,keepnear,1,0,jsubset,insave)

        pnear(pbody)=nnearbod+1
        nnearlis(pbody)=insave
        nnearbod=nnearbod+insave

        IF(nnearbod.GT.maxnearb) 
     &     CALL terror(' array overflow in reduceh ')
C               ------

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 i=1,insave
           nearbods(nnearbod-insave+i)=nearpb(isubset(jsubset(i)))
 40     CONTINUE

C   Consider neighbors within the range 2*hsmooth --> 2*hsmold.
C   -----------------------------------------------------------

        fourhsm2=four*hsmooth(pbody)*hsmooth(pbody)

        CALL WHENFGT(inear,tempvect,1,fourhsm2,isubset,ncheck)

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 50 i=1,ncheck
           jsubset(i)=nearpb(isubset(i))
           keepcold=hsmooth(jsubset(i)).GT.hsmold.OR.
     &                 ((hsmooth(jsubset(i)).EQ.hsmold).AND.
     &                 jsubset(i).LT.pbody)
           keepcrit=(templist(jsubset(i)).LT.0).AND.
     &                 (tempvect(isubset(i)).LE.four*
     &                 hsmooth(jsubset(i))**2).AND.(.NOT.keepcold)
           IF(keepcrit) THEN
              keepnear(i)=2
           ELSE
              keepnear(i)= -2
           ENDIF
 50     CONTINUE

        CALL WHENIGT(ncheck,keepnear,1,0,jsubset,insave)

        nnearlis(pbody)=nnearlis(pbody)+insave
        nnearbod=nnearbod+insave

        IF(nnearbod.GT.maxnearb) 
     &     CALL terror(' array overflow in reduceh ')
C               ------

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 60 i=1,insave
           nearbods(nnearbod-insave+i)=nearpb(isubset(jsubset(i)))
 60     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE enlargeh(pbody,nnearbod)
C
C
C***********************************************************************
C
C
C     Subroutine to increase smoothing length of body pbody so that
C     the smoothing volume contains roughly nsmooth neighbors.  The
C     neighbor list is returned from the subroutine nearpart in the
C     vector nearlist, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER i,pbody,nhtest,nearlist(nbodsmax),nnearbod,
     &          jsubset(nbodsmax),insave,npnear,keepnear(nbodsmax),
     &          nhold
        LOGICAL keepcrit,keepcold
        REAL hmin,hmax,htest,hsmold,fourhsm2,testnn,hsearch,dx,dy,dz

        EQUIVALENCE (nearlist(1),parent(1)),(jsubset(1),subindex(1)),
     &              (keepnear(1),asubp(1))

C=======================================================================

        hsearch=2.0
        npnear=0

 5      CONTINUE

        IF(npnear.LT.nsmooth) THEN

           hsearch=1.5*hsearch

           IF(hsearch.GT.1.e5)
     &        CALL terror(' search error in enlargeh ')
C                  ------

           CALL nearpart(hsearch,pbody,npnear)
C               --------

           GO TO 5

        ENDIF
        
        nhtest=npnear

        DO 10 i=1,npnear
           dx=pos(pbody,1)-pos(nearlist(i),1)
           dy=pos(pbody,2)-pos(nearlist(i),2)
           dz=pos(pbody,3)-pos(nearlist(i),3)
           IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
           IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
           IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
           IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
           IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
           IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
           tempvect(i)=dx**2+dy**2+dz**2
           isubset(i)=i
 10     CONTINUE

        testnn=ABS(REAL(npnear-nsmooth))/REAL(nsmooth)

        IF(npnear.LE.nsmooth.OR.testnn.LE.nsmtol) THEN

           htest=hsearch*hsmooth(pbody)

        ELSE

           hmin=0.
           hmax=hsearch*hsmooth(pbody)

           nhold=-1

 20        CONTINUE

           IF(ABS(REAL(nhtest-nsmooth)).GT.nsmtol*REAL(nsmooth)) THEN

              htest=0.5*(hmin+hmax)
              fourhsm2=four*htest*htest

              CALL WHENFLT(npnear,tempvect,1,fourhsm2,isubset,nhtest)

C              IF(nhtest.EQ.nhold)
C     &           CALL terror(' iteration error in enlargeh ')
C                     ------

              nhold=nhtest

              IF(nhtest.GT.nsmooth) THEN
                 hmax=htest
              ELSE
                 hmin=htest
              ENDIF

              GO TO 20

           ENDIF

        ENDIF

        nntot=nntot+nhtest
        IF(nhtest.LT.nnmin) nnmin=nhtest
        IF(nhtest.GT.nnmax) nnmax=nhtest
        nnear(pbody)=nhtest
        hsmold=hsmooth(pbody)
        hsmooth(pbody)=htest
        templist(pbody)=-1

C-----------------------------------------------------------------------
C   Select bodies to add to neighbor list.  Apply standard selection
C   criteria but do not add neighbors which have been processed and
C   have previously placed this pair on the neighbor list.
C-----------------------------------------------------------------------
        
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nhtest
           jsubset(i)=nearlist(isubset(i))
           keepcold=hsmooth(jsubset(i)).GT.hsmold.OR.
     &                 ((hsmooth(jsubset(i)).EQ.hsmold).AND.
     &                 jsubset(i).LT.pbody)
           keepcrit=((hsmooth(pbody).GT.hsmooth(jsubset(i)).OR.
     &                 ((hsmooth(pbody).EQ.hsmooth(jsubset(i))).AND.
     &                 pbody.LT.jsubset(i))).AND.templist(jsubset(i))
     &                 .GT.0).OR.((templist(jsubset(i)).LT.0).AND.
     &                 (.NOT.keepcold)).OR.((templist(jsubset(i))
     &                 .LT.0).AND.(tempvect(isubset(i)).GT.
     &                 four*hsmooth(jsubset(i))*hsmooth(jsubset(i))))
           IF(keepcrit) THEN
              keepnear(i)= 2
           ELSE
              keepnear(i)= -2
           ENDIF
 30     CONTINUE

        CALL WHENIGT(nhtest,keepnear,1,0,jsubset,insave)

        pnear(pbody)=nnearbod+1
        nnearlis(pbody)=insave
        nnearbod=nnearbod+insave

        IF(nnearbod.GT.maxnearb) 
     &     CALL terror(' array overflow in enlargeh ')
C               ------

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 i=1,insave
           nearbods(nnearbod-insave+i)=nearlist(isubset(jsubset(i)))
 40     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                SUBROUTINE nearpart(hsearch,pbody,npnear)
C
C
C***********************************************************************
C
C
C     Subroutine to search for nearest neighbors of body pbody within
C     hsearch smoothing lengths.  Vectorization is achieved by 
C     processing all cells at a given level in the tree simultaneously.
C     The neighbor list is passed back to the calling routine in the
C     vector nearlist, which is equivalenced to the common array
C     parent.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER pbody,i,nsubdiv,nearlist(nbodsmax),npnear,nnodes,nkeep,
     &          k,ibody,jsubset(nbodsmax),njsubset,nodelist(ncells),
     &          keepnear(nbodsmax),keepstak(nbodsmax)
        LOGICAL testcrit
        REAL pbottom(ndim),ptop(ndim),hsearch,sradius,dx,dy,dz,xnode,
     &       ynode,znode

        EQUIVALENCE (nearlist(1),parent(1)),(jsubset(1),subindex(1)),
     &              (nodelist(1),celllist(1)),(keepnear(1),subindex(1)),
     &              (keepstak(1),asubp(1))

C=======================================================================

        npnear=0
        sradius=hsearch*hsmooth(pbody)

        DO 10 k=1,3
           ptop(k)=pos(pbody,k)+sradius
           pbottom(k)=pos(pbody,k)-sradius
 10     CONTINUE

        nnodes=intone
        nodelist(1)=root

 20     CONTINUE

        IF(nnodes.GT.0) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nnodes
              xnode=bottom(nodelist(i),1)+0.5*cellsize(nodelist(i))
              ynode=bottom(nodelist(i),2)+0.5*cellsize(nodelist(i))
              znode=bottom(nodelist(i),3)+0.5*cellsize(nodelist(i))
              dx=0.0
              dy=0.0
              dz=0.0
              IF(xnode-pos(pbody,1).GE.hboxsize.AND.bwrap) dx=-pboxsize
              IF(xnode-pos(pbody,1).LT.-hboxsize.AND.bwrap) dx=pboxsize
              IF(ynode-pos(pbody,2).GE.hboxsize.AND.bwrap) dy=-pboxsize
              IF(ynode-pos(pbody,2).LT.-hboxsize.AND.bwrap) dy=pboxsize
              IF(znode-pos(pbody,3).GE.hboxsize.AND.bwrap) dz=-pboxsize
              IF(znode-pos(pbody,3).LT.-hboxsize.AND.bwrap) dz=pboxsize
              testcrit=     pbottom(1).LE.(bottom(nodelist(i),1)+
     &                                   cellsize(nodelist(i))+dx)
     &                 .AND.pbottom(2).LE.(bottom(nodelist(i),2)+
     &                                   cellsize(nodelist(i))+dy)
     &                 .AND.pbottom(3).LE.(bottom(nodelist(i),3)+
     &                                   cellsize(nodelist(i))+dz)
     &                 .AND.ptop(1).GE.(bottom(nodelist(i),1)+dx)
     &                 .AND.ptop(2).GE.(bottom(nodelist(i),2)+dy)
     &                 .AND.ptop(3).GE.(bottom(nodelist(i),3)+dz)
     &                 .AND.nodelist(i).NE.pbody 
              IF(testcrit.AND.nodelist(i).LE.nbodsmax) THEN
                 keepnear(i)=2
              ELSE
                 keepnear(i)= -2
              ENDIF
              IF(testcrit.AND.nodelist(i).GT.nbodsmax) THEN
                 keepstak(i)=2
              ELSE
                 keepstak(i)= -2
              ENDIF
 30        CONTINUE

           CALL WHENIGT(nnodes,keepnear,1,0,isubset,nkeep)

           DO 35 i=1,nkeep
              ibody=nodelist(isubset(i))
              dx=pos(pbody,1)-pos(ibody,1)
              dy=pos(pbody,2)-pos(ibody,2)
              dz=pos(pbody,3)-pos(ibody,3)
              IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
              IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
              IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
              IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
              IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
              IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
              tempvect(i)=dx**2+dy**2+dz**2
 35        CONTINUE

           CALL WHENFLT(nkeep,tempvect,1,sradius**2,jsubset,njsubset)

           nkeep=njsubset

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nkeep
              nearlist(npnear+i)=nodelist(isubset(jsubset(i)))
 40        CONTINUE

           npnear=npnear+nkeep

           CALL WHENIGT(nnodes,keepstak,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.npart)THEN
	      WRITE(uterm,*)8*nsubdiv
              CALL terror(' asubp overflow in nearpart ')
C                  ------
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 50        CONTINUE

           CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

           IF(nnodes.GT.nbodsmax.OR.nnodes.GT.ncells)THEN
	      WRITE(uterm,*)8*nsubdiv,nnodes
              CALL terror(' nodelist overflow in nearpart ')
C                  ------
           ENDIF

           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 20

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE celleps
C
C
C***********************************************************************
C
C
C     Subroutine to compute gravitational softening lengths of
C     cells, processing cells  in order of increasing size.  The 
C     permutation vector is stored in the common variable celllist.
C     Vectorization is achieved by simultaneously processing all 
C     cells at the same level in the hierarchy.  The softening length 
C     for each cell is a mass-weighted mean of the softening lengths 
C     of its descendents.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,fcell,lcell,i,j,nnodes,lcf1

C=======================================================================
        
C   Generate permutation of cells, according to cellsize.
C   -----------------------------------------------------

        DO 5 i=1,incells
           celllist(i)=nbodsmax+incells-(i-1)
 5      CONTINUE

C   Initialize properties of cells.
C   -------------------------------

        DO 10 p=nbodsmax+1,nbodsmax+incells
           epsvect(p)=0.0
 10     CONTINUE

C   Process cells in order of increasing size.
C   ------------------------------------------

        fcell=1

 40     CONTINUE

        IF(fcell.LE.incells) THEN

C   Determine which cells to process.
C   ---------------------------------

           DO 50 i=fcell,incells
              IF(ABS(cellsize(celllist(i))-cellsize(celllist(fcell)))
     &           .LT.0.01*cellsize(celllist(fcell))) THEN

                 lcell=i
              ELSE
                 GO TO 60
              ENDIF
 50        CONTINUE                    

 60        CONTINUE

           lcf1=lcell-fcell+1

           IF(lcf1.GT.ncells.OR.lcf1.GT.nbodsmax)
     &        CALL terror(' lcf1 overflow in celleps ')
C                  ------

C   Compute properties of the selected cells, looping over subcells.
C   ----------------------------------------------------------------

           DO 110 j=1,nsubcell

              DO 70 i=fcell,lcell
                 asubp(i-fcell+1)=subp(celllist(i),j)
 70           CONTINUE

              CALL WHENIGT(lcell-fcell+1,asubp,1,0,isubset,nnodes)

              IF(nnodes.GT.ncells.OR.nnodes.GT.nbodsmax)
     &           CALL terror(' array overflow in celleps ')
C                     ------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nnodes
                 parent(i)=celllist(isubset(i)+fcell-1)
                 asubp(i)=subp(parent(i),j)
                 epsvect(parent(i))=epsvect(parent(i))+mass(asubp(i))*
     &                              epsvect(asubp(i))
 80           CONTINUE

 110       CONTINUE


CVD$ NODEPCHK
CDIR$ IVDEP
           DO 120 i=fcell,lcell
              epsvect(celllist(i))=epsvect(celllist(i))/
     &           mass(celllist(i))
 120       CONTINUE

           fcell=lcell+1

           GO TO 40

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE accgrav(option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the gravitational acceleration for all of
C     the bodies.  Vectorization is achieved by processing all of the
C     cells at a given level in the tree simultaneously.  The local
C     variable option indicates whether the code is to compute the
C     potential and/or acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*4 option
        INTEGER p,i,j,nterms,nedge,nnonedge

C=======================================================================

C   Initialize the interaction list diagnostics.
C   --------------------------------------------
        nttot=0
        ntmin=nbodies
        ntmax=0

C   Main loop over all bodies.
C   --------------------------

        DO 100 i=1,npactive

           p=pactive(i)

C   Establish interaction lists.
C   ----------------------------
           IF(.NOT.dgravsum) THEN

              CALL treewalk(p,nterms)
C                  --------
           ELSE

              nterms=nbodies

              DO 50 j=1,nbodies
                 bodlist(j)=j
 50           CONTINUE

           ENDIF

C   Compute potential and acceleration.
C   -----------------------------------

           IF(boundary.EQ.'vacuum') THEN

              CALL gravsum(p,nterms,option)
C                  -------
           ELSE

              IF(boundary.EQ.'qperiodic') THEN

                 CALL findedge(p,nterms,nedge,nnonedge)
C                     --------
                 nterms=nnonedge+8*nedge

                 CALL qpgsum(p,nnonedge,option)
C                     ------
                 CALL edgesum(p,nedge,nnonedge,option)
C                     -------
              ELSE

                 CALL ewaldsum(p,nterms,option)
C                     --------
              ENDIF
           ENDIF                 

C   Update diagnostics, subtracting self-interaction term.
C   ------------------------------------------------------
           nterms=nterms-1
           nttot=nttot+nterms
           IF(nterms.LT.ntmin) ntmin=nterms
           IF(nterms.GT.ntmax) ntmax=nterms

 100    CONTINUE

C   Compute average number of force terms per body.
C   -----------------------------------------------
        ntavg=nttot/npactive

        RETURN
        END
C***********************************************************************
C
C
                     SUBROUTINE treewalk(p,nterms)
C
C
C***********************************************************************
C
C
C     Subroutine to walk through the tree and accumulate the list of
C     interactions for body p.  The interaction list is passed back
C     to the calling subroutine accgrav in the vector iterms, which
C     is equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nnodes,nkeep,nsubdiv,nterms,iterms(nbodsmax),
     &          nodelist(ncells),keepterm(nbodsmax)
        LOGICAL tolcrit
        REAL dxp,dyp,dzp

        EQUIVALENCE (iterms(1),bodlist(1)),(nodelist(1),celllist(1)),
     &              (keepterm(1),parent(1))

C=======================================================================

C   Initialize list of cells to examine.
C   ------------------------------------
        nterms=0
        nnodes=intone
        nodelist(1)=root
     
 10     CONTINUE

C   Loop until no cells are left to examine.
C   ----------------------------------------
        IF(nnodes.GT.0) THEN

C   Apply tolerance criterion to list of cells.
C   -------------------------------------------

           IF(boundary.EQ.'vacuum') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 20 i=1,nnodes
                 tolcrit=(tol2*((pos(p,1)-pos(nodelist(i),1))**2+
     &                   (pos(p,2)-pos(nodelist(i),2))**2+
     &                   (pos(p,3)-pos(nodelist(i),3))**2)).GE.
     &                    cellsize(nodelist(i))**2
                 IF(tolcrit) THEN
                    keepterm(i)=2
                 ELSE
                    keepterm(i) = -2
                 ENDIF
 20           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 25 i=1,nnodes
                 dxp=pos(p,1)-pos(nodelist(i),1)
                 dyp=pos(p,2)-pos(nodelist(i),2)
                 dzp=pos(p,3)-pos(nodelist(i),3)
                 IF(dxp.GE.hboxsize) dxp=dxp-pboxsize
                 IF(dxp.LT.-hboxsize) dxp=dxp+pboxsize
                 IF(dyp.GE.hboxsize) dyp=dyp-pboxsize
                 IF(dyp.LT.-hboxsize) dyp=dyp+pboxsize
                 IF(dzp.GE.hboxsize) dzp=dzp-pboxsize
                 IF(dzp.LT.-hboxsize) dzp=dzp+pboxsize
                 tolcrit=(tol2*(dxp**2+dyp**2+dzp**2)).GE.
     &                    cellsize(nodelist(i))**2
                 IF(tolcrit) THEN
                    keepterm(i)=2
                 ELSE
                    keepterm(i)= -2
                 ENDIF
 25           CONTINUE

           ENDIF

C-----------------------------------------------------------------------
C   Add cells which satisfy criterion to interaction list.  Note that,
C   depending on theta, self-interaction term may be included.
C-----------------------------------------------------------------------

           CALL WHENIGT(nnodes,keepterm,1,0,isubset,nkeep)

           IF(nterms+nkeep.GT.nbodsmax) 
     &        CALL terror(' array overflow in treewalk ')
C                  ------

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nkeep
              iterms(nterms+i)=nodelist(isubset(i))
 30        CONTINUE

           nterms=nterms+nkeep

C-----------------------------------------------------------------------
C   Add subcells of cells which fail tolerance criterion to list of
C   cells to examine.
C-----------------------------------------------------------------------

           CALL WHENILT(nnodes,keepterm,1,0,isubset,nsubdiv)

           IF(8*nsubdiv.GT.npart)THEN
	      WRITE(uterm,*)8*nsubdiv,nnodes
              CALL terror(' asubp overflow in treewalk ')
C                  ------
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nsubdiv
              asubp(i)=subp(nodelist(isubset(i)),1)
              asubp(i+nsubdiv)=subp(nodelist(isubset(i)),2)
              asubp(i+2*nsubdiv)=subp(nodelist(isubset(i)),3)
              asubp(i+3*nsubdiv)=subp(nodelist(isubset(i)),4)
              asubp(i+4*nsubdiv)=subp(nodelist(isubset(i)),5)
              asubp(i+5*nsubdiv)=subp(nodelist(isubset(i)),6)
              asubp(i+6*nsubdiv)=subp(nodelist(isubset(i)),7)
              asubp(i+7*nsubdiv)=subp(nodelist(isubset(i)),8)
 40        CONTINUE

           IF(nsubdiv.GT.0) THEN
              CALL WHENNE(8*nsubdiv,asubp,1,0,isubset,nnodes)

	      IF(nnodes.GT.nbodsmax.OR.nnodes.GT.ncells)THEN
		  WRITE(uterm,*)8*nsubdiv,nnodes
		  CALL terror(' nodelist overflow in treewalk ')
C                      ------
	      ENDIF
           ELSE
              nnodes=intzero
           ENDIF

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nnodes
              nodelist(i)=asubp(isubset(i))
 60        CONTINUE

           GO TO 10

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE gravsum(p,nterms,option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the monopole and quadrupole contributions
C     to the potential and acceleration components for body p.  The
C     interaction list is contained in the vector iterms, which is
C     equivalenced to the common array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnterm

        PARAMETER(maxnterm=nworkvec/9)

        CHARACTER*4 option
        INTEGER p,i,qindex(nbodsmax),qterms(nbodsmax),smindex(nbodsmax),
     &          nterms,iterms(nbodsmax),nqterms
        REAL r3inveff(maxnterm),rinveff(maxnterm),drdeldrg,pmass,
     &       drdotdr(maxnterm),phsm,drsm,accsm,dx(maxnterm),
     &       dy(maxnterm),dz(maxnterm),qr5inv(maxnterm),acci,
     &       phiquad(maxnterm),sdrdotdr,r2inveff(maxnterm),esoftsum

        EQUIVALENCE (iterms(1),bodlist(1)),(qindex(1),templist(1)),
     &              (qterms(1),isubset(1)),(smindex(1),parent(1)),
     &              (dx(1),workvect(1)),(dy(1),workvect(maxnterm+1)),
     &              (dz(1),workvect(2*maxnterm+1)),(r3inveff(1),
     &              workvect(3*maxnterm+1)),(rinveff(1),
     &              workvect(4*maxnterm+1)),(drdotdr(1),
     &              workvect(5*maxnterm+1)),(qr5inv(1),
     &              workvect(6*maxnterm+1)),(phiquad(1),
     &              workvect(7*maxnterm+1)),(r2inveff(1),
     &              workvect(8*maxnterm+1))

C=======================================================================

        IF(nterms.GT.maxnterm)
     &     CALL terror(' array overflow in gravsum ')
C               ------

C-----------------------------------------------------------------------
C   Compute monopole contribution; temporarily set mass of body p to
C   zero to avoid possible self-interaction contribution.
C-----------------------------------------------------------------------

        pmass=mass(p)
        mass(p)=0.

C   Loop over interaction list.
C   ---------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nterms
           dx(i)=pos(p,1)-pos(iterms(i),1)
           dy(i)=pos(p,2)-pos(iterms(i),2)
           dz(i)=pos(p,3)-pos(iterms(i),3)
           drdotdr(i)=dx(i)**2+dy(i)**2+dz(i)**2+tiny*
     &                (epsvect(p)+epsvect(iterms(i)))**2/4.
           sdrdotdr=SQRT(drdotdr(i))
           rinveff(i)=1./sdrdotdr
           r3inveff(i)=rinveff(i)/drdotdr(i)
           drdeldrg=sdrdotdr*ninterp/(epsvect(p)+epsvect(iterms(i)))
           smindex(i)=drdeldrg
           IF(ninterp.LT.smindex(i)) smindex(i)=ninterp
           IF(one.LT.drdeldrg-smindex(i)) THEN
              drsm=one
           ELSE
              drsm=drdeldrg-smindex(i)
           ENDIF
           phsm=(1.-drsm)*phsmooth(smindex(i))+
     &          drsm*phsmooth(1+smindex(i))
           accsm=(1.-drsm)*acsmooth(smindex(i))+
     &           drsm*acsmooth(1+smindex(i))
           rinveff(i)=phsm*rinveff(i)
           r3inveff(i)=accsm*r3inveff(i)
 30     CONTINUE

        IF(option.NE.'acc ') THEN

           IF(npactive.NE.nbodies) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 40 i=1,nterms
                 phi(p)=phi(p)-mass(iterms(i))*rinveff(i)
 40           CONTINUE

           ELSE

              esoftsum=0.0

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 45 i=1,nterms
                 phi(p)=phi(p)-mass(iterms(i))*rinveff(i)
                 esoftsum=esoftsum+mass(iterms(i))*(rinveff(i)-
     &                    drdotdr(i)*r3inveff(i))
 45           CONTINUE

           ENDIF

        ENDIF

        IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nterms
              acci=mass(iterms(i))*r3inveff(i)
              acc(p,1)=acc(p,1)-dx(i)*acci
              acc(p,2)=acc(p,2)-dy(i)*acci
              acc(p,3)=acc(p,3)-dz(i)*acci
 50        CONTINUE

        ENDIF

C   Reset mass of body p.
C   ---------------------
        mass(p)=pmass
 
        IF(npactive.EQ.nbodies) esofttot=esofttot+0.5*mass(p)*esoftsum

C   If required, compute quadrupole contribution.
C   ---------------------------------------------
        IF(usequad) THEN

C   Filter out bodies.
C   ------------------

           CALL WHENIGT(nterms,iterms,1,nbodsmax,qindex,nqterms)
 
C   Compute quadrupole interaction from cells.
C   ------------------------------------------
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nqterms
              qterms(i)=iterms(qindex(i))
              r2inveff(i)=rinveff(qindex(i))*rinveff(qindex(i))
              qr5inv(i)=r3inveff(qindex(i))*r2inveff(i)
              phiquad(i)=(-.5*((dx(qindex(i))**2-dz(qindex(i))**2)*
     &              quad(qterms(i),1)+(dy(qindex(i))**2-
     &              dz(qindex(i))**2)*quad(qterms(i),4))-
     &              (dx(qindex(i))*dy(qindex(i))*quad(qterms(i),2)+
     &              dx(qindex(i))*dz(qindex(i))*quad(qterms(i),3)+
     &              dy(qindex(i))*dz(qindex(i))*quad(qterms(i),5)))*
     &              qr5inv(i)
 60        CONTINUE

           IF(option.NE.'acc ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 70 i=1,nqterms
                 phi(p)=phi(p)+phiquad(i)
 70           CONTINUE

           ENDIF

           IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 80 i=1,nqterms
                 phiquad(i)=5.*phiquad(i)*r2inveff(i)
                 acc(p,1)=acc(p,1)+dx(qindex(i))*phiquad(i)+
     &                  (dx(qindex(i))*quad(qterms(i),1)+
     &                  dy(qindex(i))*quad(qterms(i),2)+
     &                  dz(qindex(i))*quad(qterms(i),3))*qr5inv(i)
                 acc(p,2)=acc(p,2)+dy(qindex(i))*phiquad(i)+
     &                  (dy(qindex(i))*quad(qterms(i),4)+
     &                  dx(qindex(i))*quad(qterms(i),2)+
     &                  dz(qindex(i))*quad(qterms(i),5))*qr5inv(i)
                 acc(p,3)=acc(p,3)+dz(qindex(i))*phiquad(i)+
     &                  (dz(qindex(i))*(-quad(qterms(i),1)-
     &                  quad(qterms(i),4))+dx(qindex(i))*
     &                  quad(qterms(i),3)+dy(qindex(i))*
     &                  quad(qterms(i),5))*qr5inv(i)
 80           CONTINUE

           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
              SUBROUTINE findedge(p,nterms,nedge,nnonedge)
C
C
C***********************************************************************
C
C
C     Subroutine to identify the cells near the edges of the
C     centered box.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER p,i,nedge,nterms,iterms(nbodsmax),cellterm(nbodsmax),
     &          ncterms,nbterms,nnofrag,cellset(nbodsmax),icell,
     &          nnonedge,jsubset(nbodsmax)
        LOGICAL testcrit
        REAL delxb,delxt,delyb,delyt,delzb,delzt

        EQUIVALENCE (iterms(1),bodlist(1)),(cellterm(1),parent(1)),
     &              (cellset(1),asubp(1)),(jsubset(1),subindex(1))
     
C=======================================================================
     
        CALL WHENIGT(nterms,iterms,1,nbodsmax,cellterm,ncterms)
     
        DO 10 i=1,ncterms
           icell=iterms(cellterm(i))
           delxb=bottom(icell,1)-pos(p,1)
           delxt=bottom(icell,1)+cellsize(icell)-pos(p,1)
           delyb=bottom(icell,2)-pos(p,2)
           delyt=bottom(icell,2)+cellsize(icell)-pos(p,2)
           delzb=bottom(icell,3)-pos(p,3)
           delzt=bottom(icell,3)+cellsize(icell)-pos(p,3)
           testcrit    = ( (delxb.LT.SIGN(onehalf,delxb)).AND.
     &                     (delxt.GE.SIGN(onehalf,delxt)) )
     &              .OR. ( (delyb.LT.SIGN(onehalf,delyb)).AND.
     &                     (delyt.GE.SIGN(onehalf,delyt)) )
     &              .OR. ( (delzb.LT.SIGN(onehalf,delzb)).AND.
     &                     (delzt.GE.SIGN(onehalf,delzt)) )
           IF(testcrit) THEN
              jsubset(i) = 2
           ELSE
              jsubset(i) = -2
           ENDIF
 10     CONTINUE
     
        CALL WHENIGT(ncterms,jsubset,1,0,cellset,nedge)
     
        DO 20 i=1,nedge
           templist(nterms-nedge+i)=iterms(cellterm(cellset(i)))
 20     CONTINUE
     
        CALL WHENILE(nterms,iterms,1,nbodsmax,isubset,nbterms)
        CALL WHENILT(ncterms,jsubset,1,0,cellset,nnofrag)
     
        DO 30 i=1,nbterms
           templist(i)=iterms(isubset(i))
 30     CONTINUE
     
        DO 40 i=1,nnofrag
           templist(nbterms+i) = iterms(cellterm(cellset(i)))
 40     CONTINUE
     
        nnonedge=nbterms+nnofrag

        DO 50 i=1,nterms
           iterms(i)=templist(i)
 50     CONTINUE
     
        RETURN
        END
C***********************************************************************
C
C
                   SUBROUTINE qpgsum(p,nterms,option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the potential and acceleration components
C     for body p, using quasi-periodic boundary conditions, from the
C     bodies and those cells not straddling the edges of the box.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnterm

        PARAMETER(maxnterm=nworkvec/9)

        CHARACTER*4 option
        INTEGER p,i,smindex(maxnterm),nterms,iterms(maxnterm)
        REAL drdeldrg,drdotdr,phsm,drsm,accsm,dx(maxnterm),dy(maxnterm),
     &       dz(maxnterm),pmass,acci,phii(maxnterm),accx(maxnterm),
     &       accy(maxnterm),accz(maxnterm),sdrdotdr,rinveff,r3inveff

        EQUIVALENCE (iterms(1),bodlist(1)),(smindex(1),parent(1)),
     &              (dx(1),workvect(1)),(dy(1),workvect(maxnterm+1)),
     &              (dz(1),workvect(2*maxnterm+1)),(phii(1),
     &              workvect(3*maxnterm+1)),(accx(1),
     &              workvect(4*maxnterm+1)),(accy(1),
     &              workvect(5*maxnterm+1)),(accz(1),
     &              workvect(6*maxnterm+1))

C=======================================================================
     
        IF(nterms.GT.maxnterm)
     &     CALL terror(' array overflow in qpgsum ')
C               ------

C-----------------------------------------------------------------------
C   Compute monopole contribution; temporarily set mass of body p to
C   zero to avoid possible self-interaction contribution.
C-----------------------------------------------------------------------
     
        pmass=mass(p)
        mass(p)=0.
     
C   Loop over interaction list.
C   ---------------------------
     
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nterms
           dx(i)=pos(p,1)-pos(iterms(i),1)
           dy(i)=pos(p,2)-pos(iterms(i),2)
           dz(i)=pos(p,3)-pos(iterms(i),3)
           IF(dx(i).GE.hboxsize) dx(i)=dx(i)-pboxsize
           IF(dx(i).LT.-hboxsize) dx(i)=dx(i)+pboxsize
           IF(dy(i).GE.hboxsize) dy(i)=dy(i)-pboxsize
           IF(dy(i).LT.-hboxsize) dy(i)=dy(i)+pboxsize
           IF(dz(i).GE.hboxsize) dz(i)=dz(i)-pboxsize
           IF(dz(i).LT.-hboxsize) dz(i)=dz(i)+pboxsize
           drdotdr=dx(i)**2+dy(i)**2+dz(i)**2+tiny*
     &             (epsvect(p)+epsvect(iterms(i)))**2/4.
           sdrdotdr=SQRT(drdotdr)
           rinveff=1./sdrdotdr
           r3inveff=rinveff/drdotdr
           drdeldrg=sdrdotdr*ninterp/(epsvect(p)+epsvect(iterms(i)))
           smindex(i)=drdeldrg
           IF(ninterp.LT.smindex(i)) smindex(i)=ninterp
           IF(one.LT.drdeldrg-smindex(i)) THEN
              drsm=one
           ELSE
              drsm=drdeldrg-smindex(i)
           ENDIF
           phsm=(1.-drsm)*phsmooth(smindex(i))+
     &          drsm*phsmooth(1+smindex(i))
           accsm=(1.-drsm)*acsmooth(smindex(i))+
     &           drsm*acsmooth(1+smindex(i))
           rinveff=phsm*rinveff
           r3inveff=accsm*r3inveff
           phii(i)= -mass(iterms(i))*rinveff
           acci= mass(iterms(i))*r3inveff
           accx(i)= -dx(i)*acci
           accy(i)= -dy(i)*acci
           accz(i)= -dz(i)*acci
 30     CONTINUE
     
C   Reset mass of body p.
C   ---------------------
        mass(p)=pmass
     
        IF(option.NE.'acc ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nterms
              phi(p)=phi(p)+phii(i)
 50        CONTINUE

        ENDIF

        IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nterms
              acc(p,1)=acc(p,1)+accx(i)
              acc(p,2)=acc(p,2)+accy(i)
              acc(p,3)=acc(p,3)+accz(i)
 60        CONTINUE

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
              SUBROUTINE edgesum(p,nedge,nnonedge,option)
C
C
C***********************************************************************
C
C
C     Subroutine to split the cells near the edges of the centered
C     box to account for the periodic boundary conditions and perform
C     force summation.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER maxnterm

        PARAMETER(maxnterm=nworkvec/9)

        CHARACTER*4 option
        INTEGER p,i,nedge,nnonedge,iterms(maxnterm),smindex(maxnterm),
     &          nterms
        REAL dx(maxnterm),dy(maxnterm),dz(maxnterm),pmass(maxnterm),
     &       delecx,delecy,delecz,xlength,ylength,zlength,delx,dely,
     &       delz,dens,xcms1,xcms2,ycms1,ycms2,zcms1,zcms2,topposx,
     &       posxbot,topposy,posybot,topposz,poszbot,acci,drdeldrg,
     &       phsm,drsm,accsm,phii(maxnterm),sdrdotdr,r3inveff,
     &       accx(maxnterm),accy(maxnterm),accz(maxnterm),rinveff,
     &       drdotdr

        EQUIVALENCE (iterms(1),bodlist(1)),(smindex(1),parent(1)),
     &              (dx(1),workvect(1)),(dy(1),workvect(maxnterm+1)),
     &              (dz(1),workvect(2*maxnterm+1)),(phii(1),
     &              workvect(3*maxnterm+1)),(accx(1),
     &              workvect(4*maxnterm+1)),(accy(1),
     &              workvect(5*maxnterm+1)),(accz(1),
     &              workvect(6*maxnterm+1)),(pmass(1),
     &              workvect(7*maxnterm+1))
     
C=======================================================================
     
        IF(8*nedge.GT.maxnterm)
     &     CALL terror(' error in edgesum--array overflow ')
C               ------

        DO 10 i=1,nedge
           templist(i)=iterms(nnonedge+i)
 10     CONTINUE

        DO 20 i=1,nedge
           iterms(i)=templist(i)
 20     CONTINUE
     
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nedge
     
C   Compute dimensions of the cloud to share.
C   -----------------------------------------
           topposx=bottom(iterms(i),1)+cellsize(iterms(i))-
     &             pos(iterms(i),1)
           posxbot=pos(iterms(i),1)-bottom(iterms(i),1)
           IF(topposx.LT.posxbot) THEN
              xlength=topposx
           ELSE
              xlength=posxbot
           ENDIF
           topposy=bottom(iterms(i),2)+cellsize(iterms(i))-
     &             pos(iterms(i),2)
           posybot=pos(iterms(i),2)-bottom(iterms(i),2)
           IF(topposy.LT.posybot) THEN
              ylength=topposy
           ELSE
              ylength=posybot
           ENDIF
           topposz=bottom(iterms(i),3)+cellsize(iterms(i))-
     &             pos(iterms(i),3)
           poszbot=pos(iterms(i),3)-bottom(iterms(i),3)
           IF(topposz.LT.poszbot) THEN
              zlength=topposz
           ELSE
              zlength=poszbot
           ENDIF
     
C   Compute displacements from the lower edge of the cloud.
C   -------------------------------------------------------
           delecx=pos(p,1)-pos(iterms(i),1)
           delecy=pos(p,2)-pos(iterms(i),2)
           delecz=pos(p,3)-pos(iterms(i),3)
           delx =( delecx+xlength-SIGN(onehalf,delecx) )*onehalf
           dely =( delecy+ylength-SIGN(onehalf,delecy) )*onehalf
           delz =( delecz+zlength-SIGN(onehalf,delecz) )*onehalf
           IF(delx.LT.zero) delx=zero
           IF(delx.GE.xlength) delx=xlength
           IF(dely.LT.zero) dely=zero
           IF(dely.GE.ylength) dely=ylength
           IF(delz.LT.zero) delz=zero
           IF(delz.GE.zlength) delz=zlength
     
C   Linearly interpolate the mass.
C   ------------------------------
           dens=mass(iterms(i))/(xlength*ylength*zlength)
           pmass(8*(i-1)+1)=delx*dely*delz*dens
           pmass(8*(i-1)+2)=delx*dely*(zlength-delz)*dens
           pmass(8*(i-1)+3)=delx*(ylength-dely)*delz*dens
           pmass(8*(i-1)+4)=delx*(ylength-dely)*(zlength-delz)*dens
           pmass(8*(i-1)+5)=(xlength-delx)*dely*delz*dens
           pmass(8*(i-1)+6)=(xlength-delx)*dely*(zlength-delz)*dens
           pmass(8*(i-1)+7)=(xlength-delx)*(ylength-dely)*delz*dens
           pmass(8*(i-1)+8)=(xlength-delx)*(ylength-dely)*(zlength-
     &                      delz)*dens
     
C   Compute the positions at which the splitted masses are affected.
C   ----------------------------------------------------------------
           xcms2=delecx-delx
           xcms1=xcms2+xlength
           IF(xcms1.LT.-hboxsize) xcms1=xcms1+pboxsize
           IF(xcms1.GE.hboxsize) xcms1=xcms1-pboxsize
           IF(xcms2.LT.-hboxsize) xcms2=xcms2+pboxsize
           IF(xcms2.GE.hboxsize) xcms2=xcms2-pboxsize
           ycms2=delecy-dely
           ycms1=ycms2+ylength
           IF(ycms1.LT.-hboxsize) ycms1=ycms1+pboxsize
           IF(ycms1.GE.hboxsize) ycms1=ycms1-pboxsize
           IF(ycms2.LT.-hboxsize) ycms2=ycms2+pboxsize
           IF(ycms2.GE.hboxsize) ycms2=ycms2-pboxsize
           zcms2=delecz-delz
           zcms1=zcms2+zlength
           IF(zcms1.LT.-hboxsize) zcms1=zcms1+pboxsize
           IF(zcms1.GE.hboxsize) zcms1=zcms1-pboxsize
           IF(zcms2.LT.-hboxsize) zcms2=zcms2+pboxsize
           IF(zcms2.GE.hboxsize) zcms2=zcms2-pboxsize
     
C   Fill the cms position arrays.
C   -----------------------------
           dx(8*(i-1)+1)=xcms1
           dy(8*(i-1)+1)=ycms1
           dz(8*(i-1)+1)=zcms1
     
           dx(8*(i-1)+2)=xcms1
           dy(8*(i-1)+2)=ycms1
           dz(8*(i-1)+2)=zcms2
     
           dx(8*(i-1)+3)=xcms1
           dy(8*(i-1)+3)=ycms2
           dz(8*(i-1)+3)=zcms1
     
           dx(8*(i-1)+4)=xcms1
           dy(8*(i-1)+4)=ycms2
           dz(8*(i-1)+4)=zcms2
     
           dx(8*(i-1)+5)=xcms2
           dy(8*(i-1)+5)=ycms1
           dz(8*(i-1)+5)=zcms1
     
           dx(8*(i-1)+6)=xcms2
           dy(8*(i-1)+6)=ycms1
           dz(8*(i-1)+6)=zcms2
     
           dx(8*(i-1)+7)=xcms2
           dy(8*(i-1)+7)=ycms2
           dz(8*(i-1)+7)=zcms1
     
           dx(8*(i-1)+8)=xcms2
           dy(8*(i-1)+8)=ycms2
           dz(8*(i-1)+8)=zcms2

 30     CONTINUE
     
        DO 40 i=1,nedge
           templist(i)=iterms(i)
 40     CONTINUE
     
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 50 i=1,nedge
           iterms(8*(i-1)+1)=templist(i)
           iterms(8*(i-1)+2)=templist(i)
           iterms(8*(i-1)+3)=templist(i)
           iterms(8*(i-1)+4)=templist(i)
           iterms(8*(i-1)+5)=templist(i)
           iterms(8*(i-1)+6)=templist(i)
           iterms(8*(i-1)+7)=templist(i)
           iterms(8*(i-1)+8)=templist(i)
 50     CONTINUE
     
        nterms=8*nedge

C   Loop over interaction list.
C   ---------------------------
     
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 60 i=1,nterms
           drdotdr=dx(i)**2+dy(i)**2+dz(i)**2+tiny*
     &             (epsvect(p)+epsvect(iterms(i)))**2/4.
           sdrdotdr=SQRT(drdotdr)
           rinveff=1./sdrdotdr
           r3inveff=rinveff/drdotdr
           drdeldrg=sdrdotdr*ninterp/(epsvect(p)+epsvect(iterms(i)))
           smindex(i)=drdeldrg
           IF(ninterp.LT.smindex(i)) smindex(i)=ninterp
           IF(one.LT.drdeldrg-smindex(i)) THEN
              drsm=one
           ELSE
              drsm=drdeldrg-smindex(i)
           ENDIF
           phsm=(1.-drsm)*phsmooth(smindex(i))+
     &          drsm*phsmooth(1+smindex(i))
           accsm=(1.-drsm)*acsmooth(smindex(i))+
     &           drsm*acsmooth(1+smindex(i))
           rinveff=phsm*rinveff
           r3inveff=accsm*r3inveff
           phii(i)= -pmass(i)*rinveff
           acci= pmass(i)*r3inveff
           accx(i)= -dx(i)*acci
           accy(i)= -dy(i)*acci
           accz(i)= -dz(i)*acci
 60     CONTINUE
     
        IF(option.NE.'acc ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 70 i=1,nterms
              phi(p)=phi(p)+phii(i)
 70        CONTINUE

        ENDIF

        IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 80 i=1,nterms
              acc(p,1)=acc(p,1)+accx(i)
              acc(p,2)=acc(p,2)+accy(i)
              acc(p,3)=acc(p,3)+accz(i)
 80        CONTINUE

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE ewaldsum(p,nterms,option)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the potential and acceleration components
C     for body p, using periodic boundary conditions.
C
C     The effect of the periodic boundary condition is taken into 
C     account as follows:
C
C        (1) The acceleration and potential correction terms are
C            tabulated on a grid using Ewald summation in subroutine
C            bcforce.  The assignment is:
C
C                x-component of acceleration correction -> bcfx
C                y-component of acceleration correction -> bcfy
C                z-component of acceleration corection  -> bcfz
C                potential correction  -> bcphi.
C
C        (2) The correction terms at the position of body p are then
C            calculated in ewaldsum by interpolating values from the
C            8 nearest grid points using cloud-in-cell assignment.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------

        INTEGER maxnterm

        PARAMETER(maxnterm=nworkvec/9)

        CHARACTER*4 option
        INTEGER p,i,smindex(maxnterm),nterms,iterms(maxnterm),nx,ny,nz,
     &          ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8
        REAL drdeldrg,drdotdr,phsm,drsm,accsm,dx(maxnterm),dy(maxnterm),
     &       dz(maxnterm),pmass,acci,xsign,ysign,zsign,x0,y0,z0,xdel,
     &       ydel,zdel,wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,fxcorr,fycorr,
     &       fzcorr,phcorr,phii(maxnterm),accx(maxnterm),accy(maxnterm),
     &       accz(maxnterm),rmaxew,sdrdotdr,rinveff,r3inveff

        EQUIVALENCE (iterms(1),bodlist(1)),(smindex(1),parent(1)),
     &              (dx(1),workvect(1)),(dy(1),workvect(maxnterm+1)),
     &              (dz(1),workvect(2*maxnterm+1)),(phii(1),
     &              workvect(3*maxnterm+1)),(accx(1),
     &              workvect(4*maxnterm+1)),(accy(1),
     &              workvect(5*maxnterm+1)),(accz(1),
     &              workvect(6*maxnterm+1))

C=======================================================================
     
        IF(nterms.GT.maxnterm)
     &     CALL terror(' array overflow in ewaldsum ')
C               ------

C-----------------------------------------------------------------------
C   Compute monopole contribution; temporarily set mass of body p to
C   zero to avoid possible self-interaction contribution.
C-----------------------------------------------------------------------
     
        pmass=mass(p)
        mass(p)=0.
     
C   Loop over interaction list.
C   ---------------------------
     
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 30 i=1,nterms
           dx(i)=pos(p,1)-pos(iterms(i),1)
           dy(i)=pos(p,2)-pos(iterms(i),2)
           dz(i)=pos(p,3)-pos(iterms(i),3)
           IF(dx(i).GE.hboxsize) dx(i)=dx(i)-pboxsize
           IF(dx(i).LT.-hboxsize) dx(i)=dx(i)+pboxsize
           IF(dy(i).GE.hboxsize) dy(i)=dy(i)-pboxsize
           IF(dy(i).LT.-hboxsize) dy(i)=dy(i)+pboxsize
           IF(dz(i).GE.hboxsize) dz(i)=dz(i)-pboxsize
           IF(dz(i).LT.-hboxsize) dz(i)=dz(i)+pboxsize
           drdotdr=dx(i)**2+dy(i)**2+dz(i)**2+tiny*
     &             (epsvect(p)+epsvect(iterms(i)))**2/4.
           sdrdotdr=SQRT(drdotdr)
           rinveff=1./sdrdotdr
           r3inveff=rinveff/drdotdr
           drdeldrg=sdrdotdr*ninterp/(epsvect(p)+epsvect(iterms(i)))
           smindex(i)=drdeldrg
           IF(ninterp.LT.smindex(i)) smindex(i)=ninterp
           IF(one.LT.drdeldrg-smindex(i)) THEN
              drsm=one
           ELSE
              drsm=drdeldrg-smindex(i)
           ENDIF
           phsm=(1.-drsm)*phsmooth(smindex(i))+
     &          drsm*phsmooth(1+smindex(i))
           accsm=(1.-drsm)*acsmooth(smindex(i))+
     &           drsm*acsmooth(1+smindex(i))
           rinveff=phsm*rinveff
           r3inveff=accsm*r3inveff
           phii(i)= -mass(iterms(i))*rinveff
           acci= mass(iterms(i))*r3inveff
           accx(i)= -dx(i)*acci
           accy(i)= -dy(i)*acci
           accz(i)= -dz(i)*acci
 30     CONTINUE
     
C   Correction terms interpolated from 8 nearest grids.
C   ---------------------------------------------------
     
CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 i=1,nterms
     
           IF(dx(i).GE.zero) THEN
              xsign=-one
           ELSE
              xsign=one
           ENDIF
           IF(dy(i).GE.zero) THEN
              ysign=-one
           ELSE
              ysign=one
           ENDIF
           IF(dz(i).GE.zero) THEN
              zsign=-one
           ELSE
              zsign=one
           ENDIF

           x0=ABS(dx(i))
           y0=ABS(dy(i))
           z0=ABS(dz(i))

           rmaxew = 2.0*REAL(nmaxew)
           nx   = INT(x0*rmaxew)
           ny   = INT(y0*rmaxew)
           nz   = INT(z0*rmaxew)
           xdel = x0 - REAL(nx)/rmaxew
           ydel = y0 - REAL(ny)/rmaxew
           zdel = z0 - REAL(nz)/rmaxew
     
           xdel = xdel * rmaxew
           ydel = ydel * rmaxew
           zdel = zdel * rmaxew

C   Weight functions.
C   -----------------
           wt1  = xdel*ydel*zdel
           wt2  = xdel*ydel*(one-zdel)
           wt3  = xdel*(one-ydel)*zdel
           wt4  = xdel*(one-ydel)*(one-zdel)
           wt5  = (one-xdel)*ydel*zdel
           wt6  = (one-xdel)*ydel*(one-zdel)
           wt7  = (one-xdel)*(one-ydel)*zdel
           wt8  = (one-xdel)*(one-ydel)*(one-zdel)
     
C   Indices of correction functions.
C   --------------------------------
           ind1=(nx+1)+nmaxew2*(ny+1)+nmaxew2s*(nz+1)
           ind2=(nx+1)+nmaxew2*(ny+1)+nmaxew2s*nz
           ind3=(nx+1)+nmaxew2*ny+nmaxew2s*(nz+1)
           ind4=(nx+1)+nmaxew2*ny+nmaxew2s*nz
           ind5=nx+nmaxew2*(ny+1)+nmaxew2s*(nz+1)
           ind6=nx+nmaxew2*(ny+1)+nmaxew2s*nz
           ind7=nx+nmaxew2*ny+nmaxew2s*(nz+1)
           ind8=nx+nmaxew2*ny+nmaxew2s*nz
     
           fxcorr=wt1*bcfxvec(ind1)+wt2*bcfxvec(ind2)+wt3*bcfxvec(ind3)+
     &            wt4*bcfxvec(ind4)+wt5*bcfxvec(ind5)+wt6*bcfxvec(ind6)+
     &            wt7*bcfxvec(ind7)+wt8*bcfxvec(ind8)
           fycorr=wt1*bcfyvec(ind1)+wt2*bcfyvec(ind2)+wt3*bcfyvec(ind3)+
     &            wt4*bcfyvec(ind4)+wt5*bcfyvec(ind5)+wt6*bcfyvec(ind6)+
     &            wt7*bcfyvec(ind7)+wt8*bcfyvec(ind8)
           fzcorr=wt1*bcfzvec(ind1)+wt2*bcfzvec(ind2)+wt3*bcfzvec(ind3)+
     &            wt4*bcfzvec(ind4)+wt5*bcfzvec(ind5)+wt6*bcfzvec(ind6)+
     &            wt7*bcfzvec(ind7)+wt8*bcfzvec(ind8)
           phcorr=wt1*bcphivec(ind1)+wt2*bcphivec(ind2)+
     &            wt3*bcphivec(ind3)+wt4*bcphivec(ind4)+
     &            wt5*bcphivec(ind5)+wt6*bcphivec(ind6)+
     &            wt7*bcphivec(ind7)+wt8*bcphivec(ind8)
     
           phii(i)  = phii(i) + mass(iterms(i))*phcorr
           accx(i)  = accx(i) + mass(iterms(i))*fxcorr*xsign
           accy(i)  = accy(i) + mass(iterms(i))*fycorr*ysign
           accz(i)  = accz(i) + mass(iterms(i))*fzcorr*zsign
     
 40     CONTINUE
     
C   Reset mass of body p.
C   ---------------------
        mass(p)=pmass

        IF(option.NE.'acc ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nterms
              phi(p)=phi(p)+phii(i)
 50        CONTINUE

        ENDIF

        IF(option.NE.'pot ') THEN

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,nterms
              acc(p,1)=acc(p,1)+accx(i)
              acc(p,2)=acc(p,2)+accy(i)
              acc(p,3)=acc(p,3)+accz(i)
 60        CONTINUE

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                SUBROUTINE srchbox(p,npnear,pbottom,ptop)
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the search box for the nearest
C     neighbor detection for a group of particles specified by
C     the argument p.  The list of neighbors is returned through
C     the vector nearlist, which is equivalenced to the common
C     array bodlist.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,j,k,nearlist(nbodsmax),npgroup,ismin,ismax,npnear
        REAL pbottom(ndim),ptop(ndim),rsearch

        EQUIVALENCE (nearlist(1),bodlist(1))

C=======================================================================
 
        npgroup=pgroupb(p+1)-pgroupb(p)

        DO 30 k=1,3

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 10 j=pgroupb(p),pgroupb(p+1)-1
              IF(groupbod(j).LE.nsph) THEN
                 rsearch=hsmooth(groupbod(j))
              ELSE
                 rsearch=epsvect(groupbod(j))
              ENDIF
              tempvect(j-pgroupb(p)+1)=pos(groupbod(j),k)-2.*rsearch
 10        CONTINUE

           pbottom(k)=tempvect(ISMIN(npgroup,tempvect,1))

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 20 j=pgroupb(p),pgroupb(p+1)-1
              IF(groupbod(j).LE.nsph) THEN
                 rsearch=hsmooth(groupbod(j))
              ELSE
                 rsearch=epsvect(groupbod(j))
              ENDIF
              tempvect(j-pgroupb(p)+1)=pos(groupbod(j),k)+2.*rsearch
 20        CONTINUE

           ptop(k)=tempvect(ISMAX(npgroup,tempvect,1))

 30     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 40 j=pgroupb(p),pgroupb(p+1)-1
           nearlist(npnear+j-pgroupb(p)+1)=groupbod(j)
 40     CONTINUE

        npnear=npnear+npgroup

        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outstate(n)
C
C
C***********************************************************************
C
C
C     Subroutine to output information about the system state to
C     the log and body data files.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER n,ioutcosm

        SAVE ioutcosm

        DATA ioutcosm/1/

C=======================================================================

        CALL outterm(' step completed: ',n)
C            -------

        IF(n.EQ.0) THEN

           CALL outlog(0)
C               ------
           CALL energy
C               ------
           CALL outbods
C               -------

        ELSE

           IF(MOD(n,noutlog).EQ.0.OR.(MOD(n,noutbod).EQ.0.AND.(.NOT.
     &        comove)).OR.(comove.AND.expanpar.GE.exparout(ioutcosm)))
     &     THEN

              IF((MOD(n,noutbod).EQ.0.AND.(.NOT.comove)).OR.(comove.AND.
     &           expanpar.GE.exparout(ioutcosm))) CALL outdump(0)
C                                                      -------

              IF(comove.AND.exparout(ioutcosm).LE.zero) THEN
                 CALL endrun
C                     ------
                 STOP
              ENDIF

              CALL corrpos('correct')
C                  -------
              CALL zeropot
C                  -------
              CALL gravity('pot ')
C                  -------

              IF(MOD(n,noutlog).EQ.0) THEN

                 CALL outlog(n)
C                     ------
                 CALL energy
C                     ------
              ENDIF

              IF((MOD(n,noutbod).EQ.0.AND.(.NOT.comove)).OR.(comove.AND.
     &           expanpar.GE.exparout(ioutcosm)))
     &        THEN 
                 CALL outbods
C                     -------
                 IF(comove) THEN
                    ioutcosm=ioutcosm+1

                    IF(ioutcosm.GT.noutcosm) THEN
                       CALL endrun
C                           ------
                       STOP
                    ENDIF

                    IF(ioutcosm.GT.maxnoutc) 
     &                 CALL terror(' i/o error in outstate ')
C                           ------
                 ENDIF
              ENDIF

              CALL corrpos('reset  ')
C                  -------
           ENDIF

        ENDIF
 
        RETURN
        END
C***********************************************************************
C
C
                        SUBROUTINE outlog(istep)
C
C
C***********************************************************************
C
C
C     Subroutine to monitor status of the program by writing to
C     the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER istep,i,tstepbin(20),p
        REAL tistep

C=======================================================================
 
C   If first call, write header.
C   ----------------------------
        IF(istep.EQ.0) CALL outhead
C                           -------
 
C-----------------------------------------------------------------------
C   Output system time and force evaluation diagnostics.
C-----------------------------------------------------------------------

        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)
 
        IF(cosmo) THEN
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,70,REC=ireclog) tnow,expanpar,redshift,incellsg
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
        ELSE
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,75,REC=ireclog) tnow,incellsg
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
        ENDIF

        ireclog=ireclog+1
        WRITE(ulog,80,REC=ireclog) nttot,ntmin,ntmax,ntavg
        ireclog=ireclog+1
        WRITE(ulog,90,REC=ireclog) nntot,nnmin,nnmax,nnavg

        IF(variabls) THEN
           ireclog=ireclog+1
           WRITE(ulog,95,REC=ireclog) nstot,nsmin,nsmax,nsavg
        ENDIF

        ireclog=ireclog+1 
        WRITE(ulog,65,REC=ireclog)

        IF(usesph) THEN

           DO 10 i=1,20
              tstepbin(i)=0
 10        CONTINUE

           DO 20 p=1,nsph
              tistep=itimestp(p)
              tistep=0.5+LOG(tistep)/log2
              templist(p)=INT(tistep)+1
 20        CONTINUE

           DO 30 p=1,nsph
              tstepbin(templist(p))=tstepbin(templist(p))+1
 30        CONTINUE

           ireclog=ireclog+1
           WRITE(ulog,100,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,110,REC=ireclog) (tstepbin(i),i=1,8)
           ireclog=ireclog+1
           WRITE(ulog,110,REC=ireclog) (tstepbin(i),i=9,16)

        ENDIF

        IF(nsph.NE.nbodies) THEN

           DO 40 i=1,20
              tstepbin(i)=0
 40        CONTINUE

           DO 50 p=nsph+1,nbodies
              tistep=itimestp(p)
              tistep=0.5+LOG(tistep)/log2
              templist(p)=INT(tistep)+1
 50        CONTINUE

           DO 60 p=nsph+1,nbodies
              tstepbin(templist(p))=tstepbin(templist(p))+1
 60        CONTINUE

           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,120,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,65,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,110,REC=ireclog) (tstepbin(i),i=1,8)
           ireclog=ireclog+1
           WRITE(ulog,110,REC=ireclog) (tstepbin(i),i=9,16)

        ENDIF

        CLOSE(UNIT=ulog)

 65     FORMAT(' ')
 70     FORMAT(2x,'time: ',1pe11.3,5x,'aexp: ',1pe11.3,5x,'Z: ',
     &         0pf5.2,5x,'ncells: ',1i5)
 75     FORMAT(2x,'time: ',1pe12.4,5x,'ncells: ',1i5)
 80     FORMAT(7x,'nttot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)
 90     FORMAT(7x,'nntot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)
 95     FORMAT(7x,'nstot, min, max, avg = ',1i8,5x,1i5,5x,1i5,5x,1i5)
 100    FORMAT(7x,'distribution of time steps for SPH particles : ')
 110    FORMAT(12x,8i6)
 120    FORMAT(7x,'distribution of time steps for non-SPH particles : ')

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE outhead
C
C
C***********************************************************************
C
C
C     Subroutine to output a standard header to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,25,REC=ireclog) headline
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,28,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,29,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,30,REC=ireclog) nbodies,nsteps,noutbod,noutlog
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog) dtime,eps,usequad,tol
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,42,REC=ireclog) dgravsum,variabls,nsvolume
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,45,REC=ireclog) fixthalo,nsat,selfgrav
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,47,REC=ireclog) vhalo,gamhalo,rhalo

        IF(usebh) THEN
           ireclog=ireclog+1
           WRITE(ulog,20,REC=ireclog)
           ireclog=ireclog+1
           WRITE(ulog,48,REC=ireclog) usebh,massbh,etabh,tstartbh
        ENDIF

        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,50,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,55,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,57,REC=ireclog) usesph,sphinit,variablh,variablg,
     &                             starform
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,58,REC=ireclog) nsph,courant,nsmooth,artfvisc
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,59,REC=ireclog) uentropy,consthsm,cfrict,friclucy
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,60,REC=ireclog) gamma,alpha,beta,epssph
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,61,REC=ireclog) epsgas,friction,ethinit,entinit
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,62,REC=ireclog) readinas,outputas,readinet,outputet
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,63,REC=ireclog) symmetry,geogradp,sphfiltr,epsfiltr
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,70,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,80,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,90,REC=ireclog) inittbin,mintstep,etol,ntvector
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,100,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,110,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,120,REC=ireclog) radiate,aheatmh,bheatmh2
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,130,REC=ireclog) ccoolmh2,meanmwt,mhboltz
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,140,REC=ireclog) comptmh,slowcool,ctcutoff
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,145,REC=ireclog) mintemp,fhydrogn
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,150,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,160,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,170,REC=ireclog) cosmo,comove,hubble0,omega0
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,180,REC=ireclog) boundary,noutcosm
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,190,REC=ireclog) pboxsize,ecosm
        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
 
        CLOSE(UNIT=ulog)

 10     FORMAT(1x,72('*'))
 20     FORMAT(1x,'*',70(' '),'*')
 25     FORMAT(1x,'*',10x,1a50,10x,'*')
 28     FORMAT(1x,'*',4x,'Input parameters:',49x,'*')
 29     FORMAT(1x,'*',4x,'----------------',50x,'*')
 30     FORMAT(1x,'*',8x,'nbodies=',1i7,2x,'nsteps=',1i4,4x,
     &         'noutbods=',1i4,2x,'noutlog=',1i4,3x,'*')
 40     FORMAT(1x,'*',8x,'dtime=',1pe10.3,1x,'eps=',1pe10.3,1x,
     &         'usequad=',1l1,6x,'tol=',0pf5.2,6x,'*')
 42     FORMAT(1x,'*',8x,'dgravsum=',1l1,2x,'variabls=',1l1,2x,
     &         'nsvolume=',1i7,22x,'*')
 45     FORMAT(1x,'*',8x,'fixthalo=',1l1,2x,'nsat=',1i7,2x,
     &         'selfgrav=',1l1,26x,'*')
 47     FORMAT(1x,'*',8x,'vhalo=',1pe10.3,1x,'gamhalo=',1pe10.3,
     &         1x,'rhalo=',1pe10.3,10x,'*')
 48     FORMAT(1x,'*',8x,'usebh=',1l1,'massbh=',1pe10.3,1x,
     &         'etabh=',1pe10.3,1x,'tstartbh=',1pe10.3,1x,'*')
 50     FORMAT(1x,'*',4x,'SPH parameters:',51x,'*')
 55     FORMAT(1x,'*',4x,'--------------',52x,'*')
 57     FORMAT(1x,'*',8x,'usesph=',1l1,2x,'sphinit=',1l1,2x,
     &         'variablh=',1l1,2x,'variablg=',1l1,2x,
     &         'starform=',1l1,7x,'*')
 58     FORMAT(1x,'*',8x,'nsph=',1i7,2x,'courant=',1f5.3,2x,
     &         'nsmooth=',1i4,2x,'artfvisc= ',1a4,5x,'*')
 59     FORMAT(1x,'*',8x,'uentropy=',1l1,1x,'consthsm=',1pe10.3,
     &         1x,'cfrict=',1pe10.3,1x,'friclucy=',1l1,3x,'*')
 60     FORMAT(1x,'*',8x,'gamma=',1f5.3,2x,'alpha=',1f5.3,2x,
     &         'beta=',1f5.3,2x,'epssph=',1f5.3,12x,'*')
 61     FORMAT(1x,'*',8x,'epsgas=',1pe10.3,2x,'friction=',1l1,2x,
     &         'ethinit=',1l1,2x,'entinit=',1l1,11x,'*')
 62     FORMAT(1x,'*',8x,'readinas=',1l1,2x,'outputas=',1l1,2x,
     &         'readinet=',1l1,2x,'outputet=',1l1,16x,'*')
 63     FORMAT(1x,'*',8x,'symmetry=',1a2,2x,'geogradp=',1l1,
     &         2x,'sphfiltr=',1l1,2x,'epsfiltr=',1pe10.3,6x,'*')
 70     FORMAT(1x,'*',4x,'Time step parameters:',45x,'*')
 80     FORMAT(1x,'*',4x,'--------------------',46x,'*')
 90     FORMAT(1x,'*',8x,'inittbin=',1i5,2x,'mintstep=',1i5,2x,
     &         'etol=',1f7.3,2x,'ntvector=',1i5,2x,'*')
 100    FORMAT(1x,'*',4x,'Radation parameters:',46x,'*')
 110    FORMAT(1x,'*',4x,'-------------------',47x,'*')
 120    FORMAT(1x,'*',8x,'radiate=',1l1,2x,'aheatmh=',1pe10.3,2x,
     &         'bheatmh2=',1pe10.3,12x,'*')
 130    FORMAT(1x,'*',8x,'ccoolmh2=',1pe10.3,2x,'meanmwt=',
     &         1pe10.3,2x,'mhboltz=',1pe10.3,3x,'*')
 140    FORMAT(1x,'*',8x,'comptmh=',1pe10.3,2x,'slowcool=',
     &         1pe10.3,2x,'ctcutoff=',1pe10.3,2x,'*')
 145    FORMAT(1x,'*',8x,'mintemp=',1pe10.3,2x,'fhydrogn=',
     &         1pe10.3,23x,'*')
 150    FORMAT(1x,'*',4x,'Cosmology parameters:',45x,'*')
 160    FORMAT(1x,'*',4x,'--------------------',46x,'*')
 170    FORMAT(1x,'*',8x,'cosmo=',1l1,2x,'comove=',1l1,2x,
     &         'hubble0=',1pe10.3,6x,'omega0=',0pf5.2,7x,'*')
 180    FORMAT(1x,'*',8x,'boundary= ',1a9,2x,
     &         2x,'noutcosm=',1i6,4x,'*')
 190    FORMAT(1x,'*',8x,'pboxsize=',1pe11.3,2x,'ecosm=',
     &         1pe11.3,23x,'*')

        RETURN
        END
C***********************************************************************
C
C
                            SUBROUTINE energy
C
C
C***********************************************************************
C
C
C     Subroutine to compute diagnostics for the system: total energy,
C     total kinetic energy, total potential energy, angular momentum,
C     center of mass coordinates, and center of mass velocity.  The
C     local variable p is a pointer to the bodies.  The local
C     variables cmpos, cmvel, and amvec are center of mass position 
C     and velocity, and total angular momentum of the system, 
C     respectively.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,k
        LOGICAL firstc
        REAL cmpos(ndim),cmvel(ndim),amvec(3),epold,tepold,hubold,
     &       decosm,etold,avgrold,avgpold,avgp,avgrho

        SAVE firstc,epold,tepold,hubold,etold,avgrold,avgpold,avgp,
     &       avgrho

        DATA firstc/.TRUE./

C=======================================================================

        IF(firstc) THEN
           IF(cosmo.AND.comove) THEN
              epold=zero
              tepold=zero
              hubold=hubble
C##### LH 050990
              etold=zero
              avgrold=zero
              avgpold=zero
C##### LH 050990
           ENDIF
        ENDIF

C   Zero the accumulators for system diagnostics.
C   ---------------------------------------------
        mgas=0.
        mstar=0.
        mtot=0.
        ektot=0.
        eptot=0.
        ethtot=0.

        DO 100 k=1,ndim
           cmpos(k)=0.
           cmvel(k)=0.
 100    CONTINUE

        DO 120 k=1,3
           amvec(k)=0.
 120    CONTINUE

C-----------------------------------------------------------------------
C   Loop over bodies to compute system mass and potential energy.
C-----------------------------------------------------------------------

        DO 135 p=1,nsph
           mgas=mgas+mass(p)
 135    CONTINUE

        DO 140 p=nbodies-nstar+1,nbodies
           mstar=mstar+mass(p)
 140    CONTINUE

        IF(boundary.EQ.'vacuum'.OR.boundary.EQ.'qperiodic') THEN

           DO 150 p=1,nbodies
              mtot=mtot+mass(p)
              eptot=eptot+.5*mass(p)*phi(p)+mass(p)*phiext(p)
 150       CONTINUE

        ELSE

           DO 160 p=1,nbodies
              mtot=mtot+mass(p)
              eptot=eptot+.5*mass(p)*(phi(p)+two*alphaew*mass(p)/
     &              sqrtpi)+mass(p)*phiext(p)
 160       CONTINUE

C   Include correction due to neutralizing background.
C   --------------------------------------------------

           eptot=eptot+2.0*(mtot*sqrtpi/(2.0*alphaew))**2

           IF(.NOT.selfgrav) eptot=0.0

        ENDIF

C   Compute thermal energy.
C   -----------------------

        IF(.NOT.isotherm) THEN

           DO 170 p=1,nsph
              ethtot=ethtot+mass(p)*csound(p)**2/(gamma*gamma1)
 170       CONTINUE

        ELSE

           DO 180 p=1,nsph
              ethtot=ethtot+mass(p)*csound(p)**2
 180       CONTINUE

        ENDIF

C-----------------------------------------------------------------------
C   Compute system kinetic energy, components of center of mass
C   position and velocity.
C-----------------------------------------------------------------------

        DO 250 k=1,ndim
           DO 200 p=1,nbodies
              ektot=ektot+.5*mass(p)*vel(p,k)*vel(p,k)
              cmpos(k)=cmpos(k)+mass(p)*pos(p,k)
              cmvel(k)=cmvel(k)+mass(p)*vel(p,k)
 200       CONTINUE
           cmvel(k)=cmvel(k)/mtot
           cmpos(k)=cmpos(k)/mtot
 250    CONTINUE

        IF(cosmo.AND.comove) THEN
           ektot=ektot*cosmofac**4
           eptot=cosmofac*eptot

C##### LH 050990

           ethtot=cosmofac**2*ethtot

           avgrho=zero
           avgp=zero

           DO 255 p=1,nsph
              avgrho=avgrho+hsmooth(p)**3*rho(p)
              avgp=avgp+hsmooth(p)**3*rho(p)*csound(p)**2/gamma
 255       CONTINUE

           avgp=cosmofac**2*avgp

C##### LH 050990

        ENDIF

        IF((.NOT.firstc).AND.cosmo.AND.comove) THEN
C##### LH 050990
           IF(nsph.GT.0) THEN
              decosm=0.5*((epold*hubold+eptot*hubble)+2.*(etold*hubold+
     &               ethtot*hubble)-3.*(avgp*hubble/avgrho+avgpold*
     &               hubold/avgrold))*(tnow-tepold)
           ELSE
              decosm=0.5*(epold*hubold+eptot*hubble)*(tnow-tepold)
           ENDIF
C##### LH 050990
           ecosm=ecosm+decosm
        ENDIF

C   Compute total system energy.
C   ----------------------------
        etot=ektot+eptot+ethtot

        IF(cosmo.AND.comove) etot=etot-ecosm

        IF(usebh) THEN
           DO 270 k=1,ndim
              etot=etot+mtot*cmvel(k)*velbh(k)
 270       CONTINUE
        ENDIF

C-----------------------------------------------------------------------
C   Compute angular momentum of the system.
C-----------------------------------------------------------------------

        IF(ndim.EQ.2) THEN

           DO 300 p=1,nbodies
              amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-
     &                 pos(p,2)*vel(p,1))
 300       CONTINUE

           IF(usebh) THEN
              amvecbh(3)=amvec(3)
              amvec(3)=amvec(3)+mtot*(posbh(1)*cmvel(2)-posbh(2)*
     &                 cmvel(1)+cmpos(1)*velbh(2)-cmpos(2)*velbh(1))
           ENDIF

        ELSE IF(ndim.EQ.3) THEN

           DO 400 p=1,nbodies
              amvec(1)=amvec(1)+mass(p)*(pos(p,2)*vel(p,3)-
     &                 pos(p,3)*vel(p,2))
              amvec(2)=amvec(2)+mass(p)*(pos(p,3)*vel(p,1)-
     &                 pos(p,1)*vel(p,3))
              amvec(3)=amvec(3)+mass(p)*(pos(p,1)*vel(p,2)-
     &                 pos(p,2)*vel(p,1))
 400       CONTINUE

           IF(usebh) THEN
              amvecbh(1)=amvec(1)
              amvecbh(2)=amvec(2)
              amvecbh(3)=amvec(3)
              amvec(1)=amvec(1)+mtot*(posbh(2)*cmvel(3)-posbh(3)*
     &                 cmvel(2)+cmpos(2)*velbh(3)-cmpos(3)*velbh(2))
              amvec(2)=amvec(2)+mtot*(posbh(3)*cmvel(1)-posbh(1)*
     &                 cmvel(3)+cmpos(3)*velbh(1)-cmpos(1)*velbh(3))
              amvec(3)=amvec(3)+mtot*(posbh(1)*cmvel(2)-posbh(2)*
     &                 cmvel(1)+cmpos(1)*velbh(2)-cmpos(2)*velbh(1))
           ENDIF

        ENDIF

        IF(usebh) THEN
           DO 500 k=1,ndim
              cmvelbh(k)=cmvel(k)
              cmposbh(k)=cmpos(k)
              cmvel(k)=cmvel(k)+velbh(k)
              cmpos(k)=cmpos(k)+posbh(k)
 500       CONTINUE
        ENDIF

C   Write diagnostics to the log file.
C   ----------------------------------
        CALL outenrgy(amvec,cmpos,cmvel)
C            --------

        IF(cosmo.AND.comove) THEN
           epold=eptot
           tepold=tnow
           hubold=hubble
C##### LH 050990
           etold=ethtot
           avgrold=avgrho
           avgpold=avgp
C##### LH 050990
        ENDIF

        IF(firstc) firstc=.FALSE.

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE outbods
C
C
C***********************************************************************
C
C
C     Subroutine to output the body data.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*3 sstring
        CHARACTER*7 filename
        CHARACTER*8 filepar
        CHARACTER*10 nstring
        INTEGER p,ndimo,istring,nsnap,k

        SAVE nsnap,nstring

        DATA nsnap/0/,nstring/'0123456789'/

C=======================================================================
 
        nsnap=nsnap+1

        sstring(1:1)=nstring(1+nsnap/100:1+nsnap/100)
        istring=1+MOD(nsnap,100)/10
        sstring(2:2)=nstring(istring:istring)
        istring=1+MOD(nsnap,10)
        sstring(3:3)=nstring(istring:istring)
        filepar=obodfile
        filename=filepar(1:4)//sstring(1:3)

        OPEN(UNIT=ubodsout,FILE=filename,STATUS='NEW')

        ndimo=ndim

        WRITE(ubodsout,200) nbodies,nsph,nstar
        WRITE(ubodsout,200) ndimo
        WRITE(ubodsout,210) tnow
 
        DO 10 p=1,nbodies
           WRITE(ubodsout,210) mass(p)
 10     CONTINUE

        DO 25 k=1,ndim
           DO 20 p=1,nbodies
              WRITE(ubodsout,210) pos(p,k)
 20        CONTINUE
 25     CONTINUE

        DO 35 k=1,ndim
           DO 30 p=1,nbodies
              WRITE(ubodsout,210) vel(p,k)
 30        CONTINUE
 35     CONTINUE

        DO 40 p=nsph+1,nbodies
           WRITE(ubodsout,210) epsvect(p)
 40     CONTINUE

        DO 110 p=1,nsph
           WRITE(ubodsout,210) rho(p)
 110    CONTINUE

        IF(outputet) THEN
           IF(.NOT.isotherm) THEN
              DO 50 p=1,nsph
                 tempvect(p)=csound(p)**2/(gamma*gamma1)
 50           CONTINUE
           ELSE
              DO 55 p=1,nsph
                 tempvect(p)=csound(p)**2
 55           CONTINUE
           ENDIF
        ELSE
           IF(outputas) THEN
              DO 60 p=1,nsph
                 tempvect(p)=csound(p)**2/(gamma*rho(p)**gamma1)
 60           CONTINUE
           ELSE
              DO 70 p=1,nsph
                 tempvect(p)=csound(p)**2*meanmwt*mhboltz/gamma
 70           CONTINUE
           ENDIF
        ENDIF

        DO 120 p=1,nsph
           WRITE(ubodsout,210) tempvect(p)
 120    CONTINUE

        DO 130 p=1,nsph
           WRITE(ubodsout,210) hsmooth(p)
 130    CONTINUE

        DO 140 p=1,nsph
           WRITE(ubodsout,210) metals(p)
 140    CONTINUE

        DO 150 p=nbodies-nstar+1,nbodies
           WRITE(ubodsout,210) metals(p)
 150    CONTINUE

        DO 155 p=nbodies-nstar+1,nbodies
           WRITE(ubodsout,210) tform(p)
 155    CONTINUE

        DO 160 p=1,nbodies
           tempvect(p)=phi(p)+phiext(p)
           WRITE(ubodsout,210) tempvect(p)
 160    CONTINUE

 200    FORMAT(1x,5(1i6))
 210    FORMAT(1x,10(1pe14.6))

        CLOSE(UNIT=ubodsout)

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE outdump(iopt)
C
C
C***********************************************************************
C
C
C     Subroutine to dump state of system to an ascii data file.  The
C     argument iopt indicates whether the data is output to a normal
C     SYSDUMP file (iopt=0) or a crash dump file (iopt=1).
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,iopt,udumpout

C=======================================================================
 
        IF(iopt.EQ.0) THEN
           OPEN(UNIT=uboddump,FILE=dumpfile,STATUS='UNKNOWN')
           udumpout=uboddump
        ELSE
           OPEN(UNIT=ucrash,FILE=crashfil,STATUS='NEW')
           udumpout=ucrash
        ENDIF

C   Output system state.
C   --------------------

        WRITE(udumpout,500) nbodies,nsph,nstar

        WRITE(udumpout,510) tnow,tpos
        WRITE(udumpout,*) dtime
 
        WRITE(udumpout,*) mtot
        WRITE(udumpout,*) ektot
        WRITE(udumpout,*) eptot
        WRITE(udumpout,*) snheat
        WRITE(udumpout,*) eradiate
        WRITE(udumpout,*) esofttot
        WRITE(udumpout,*) enerror

        WRITE(udumpout,*) teth
        WRITE(udumpout,*) upbin
        WRITE(udumpout,*) stime
        WRITE(udumpout,*) tsteppos
        WRITE(udumpout,*) endstep
        WRITE(udumpout,*) trad

        DO 70 p=1,nbodies
           WRITE(udumpout,510) mass(p)
 70     CONTINUE

        DO 80 p=1,nbodies
           WRITE(udumpout,510) pos(p,1),pos(p,2),pos(p,3)
 80     CONTINUE

        DO 90 p=1,nbodies
           WRITE(udumpout,510) vel(p,1),vel(p,2),vel(p,3)
 90     CONTINUE

        DO 105 p=1,nbodies
           WRITE(udumpout,510) epsvect(p)
 105    CONTINUE

        DO 100 p=1,nsph
           WRITE(udumpout,510) rho(p)
 100    CONTINUE

        DO 110 p=1,nsph
           WRITE(udumpout,510) ethermal(p)
 110    CONTINUE

        DO 120 p=1,nsph
           WRITE(udumpout,510) hsmooth(p)
 120    CONTINUE

        DO 130 p=1,nsph
           WRITE(udumpout,510) csound(p)
 130    CONTINUE

        DO 140 p=1,nsph
           WRITE(udumpout,510) ethold(p)
 140    CONTINUE

        DO 150 p=1,nsph
           WRITE(udumpout,510) dethdt(p)
 150    CONTINUE

        DO 160 p=1,nsph
           WRITE(udumpout,510) dethold(p)
 160    CONTINUE

        DO 170 p=1,nsph
           WRITE(udumpout,510) tvel(p)
 170    CONTINUE

        DO 180 p=1,nsph
           WRITE(udumpout,510) mumaxdvh(p)
 180    CONTINUE

        DO 190 p=1,nsph
           WRITE(udumpout,510) hsmdivv(p)
 190    CONTINUE

        DO 200 p=1,nsph
           WRITE(udumpout,510) derad(p)
 200    CONTINUE

        DO 204 p=1,nsph
           WRITE(udumpout,510) hsmcurlv(p)
 204    CONTINUE

        DO 210 p=1,nbodies
           WRITE(udumpout,500) nnear(p)
 210    CONTINUE

        DO 220 p=1,nbodies
           WRITE(udumpout,510) phi(p),phiext(p)
 220    CONTINUE

        DO 230 p=1,nbodies
           WRITE(udumpout,510) acc(p,1),acc(p,2),acc(p,3)
 230    CONTINUE

        DO 240 p=1,nbodies
           WRITE(udumpout,500) itimestp(p)
 240    CONTINUE

        DO 250 p=1,nbodies
           WRITE(udumpout,500) otimestp(p)
 250    CONTINUE

        DO 260 p=1,nsph
           WRITE(udumpout,510) metals(p)
 260    CONTINUE

        DO 270 p=nbodies-nstar+1,nbodies
           WRITE(udumpout,510) metals(p)
 270    CONTINUE

        DO 280 p=nbodies-nstar+1,nbodies
           WRITE(udumpout,510) tform(p)
 280    CONTINUE

        CLOSE(UNIT=udumpout)

 500    FORMAT(1x,10(1i10))
 510    FORMAT(1x,10(1pe22.14))

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE corrpos(rc)
C
C
C***********************************************************************
C
C
C     Subroutine to apply a correction factor to the positions to
C     maintain second order accuracy when outputting particle data
C     to body data file or when computing energy diagnostics.  The
C     argument rc indicates whether the correction factor is to be
C     applied (correct) or removed (reset).
C  
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 rc
        INTEGER p,k
        REAL dt2,rcsign,acceff,dt

C=======================================================================

C   Loop over all spatial coordinates for all bodies.
C   -------------------------------------------------
  
        IF(rc.EQ.'correct') THEN
           rcsign=-1.
        ELSE
           rcsign=1.
        ENDIF

        IF(boundary.EQ.'vacuum') THEN

           IF((.NOT.friction).OR.friclucy) THEN

              DO 200 k=1,ndim
                 DO 100 p=1,nbodies
                    dt2=(dtime/itimestp(p))**2
                    pos(p,k)=pos(p,k)+rcsign*acc(p,k)*dt2/8.
 100             CONTINUE
 200          CONTINUE

           ELSE

              DO 250 k=1,ndim
                 DO 230 p=1,nbodies
                    dt=dtime/itimestp(p)
                    dt2=dt**2
                    IF(p.LE.nsph) THEN
                       acceff=cfrict/dt
                    ELSE
                       acceff=zero
                    ENDIF
                    pos(p,k)=pos(p,k)+rcsign*(acc(p,k)-acceff*
     &                       vel(p,k))*dt2/8.
 230             CONTINUE
 250          CONTINUE

           ENDIF

        ELSE

           IF(.NOT.comove) THEN

              DO 400 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 300 p=1,nbodies
                    dt2=(dtime/itimestp(p))**2
                    pos(p,k)=pos(p,k)+rcsign*acc(p,k)*dt2/8.
                    IF(pos(p,k).GE.hboxsize) pos(p,k)=pos(p,k)-pboxsize
                    IF(pos(p,k).LT.-hboxsize) pos(p,k)=pos(p,k)+pboxsize
 300             CONTINUE

 400          CONTINUE

           ELSE

              DO 600 k=1,ndim

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 500 p=1,nbodies
                    dt2=(dtime/itimestp(p))**2
                    pos(p,k)=pos(p,k)+rcsign*(acc(p,k)*dt2/
     &                       (8.*cosmof3)-vel(p,k)*dt2*hubble/4.)
                    IF(pos(p,k).GE.hboxsize) pos(p,k)=pos(p,k)-pboxsize
                    IF(pos(p,k).LT.-hboxsize) pos(p,k)=pos(p,k)+pboxsize
 500             CONTINUE
 600          CONTINUE

           ENDIF

        ENDIF

        RETURN
        END

C***********************************************************************
C
C
                           SUBROUTINE endrun
C
C
C***********************************************************************
C
C
C     Subroutine to end the simulation.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        REAL second

C=======================================================================

        cputime1=SECOND()

        CALL outcpu
C            ------
        CALL stopout
C            -------

        RETURN
        END
C***********************************************************************
C
C
                         SUBROUTINE outcpu
C
C
C***********************************************************************
C
C
C     Subroutine to output cpu timing data to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

C   Output timing data to the log file.
C   -----------------------------------
        cputime=cputime1-cputime0

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog) cputime
 
        CLOSE(UNIT=ulog)

 10     FORMAT(' ')
 20     FORMAT(' Total cpu time used (seconds) : ',1pe12.4)
 
        RETURN
        END
C***********************************************************************
C
C
                           SUBROUTINE stopout
C
C
C***********************************************************************
C
C
C     Subroutine to close the open output files.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C=======================================================================
 
C   Close the open files.
C   ---------------------
C       CLOSE(UNIT=ulog)
 
        RETURN
        END
C***********************************************************************
C
C
                         FUNCTION ismax(n,x,inc)
C
C
C***********************************************************************
C
C
C     Function to locate index of maximum element of a real vector.
C
C
C=======================================================================

        INTEGER ismax,n,inc,i
        REAL x(1),xmax

        ismax=1
        xmax=x(1)

        DO 10 i=2,n,inc
           IF(x(i).GT.xmax) THEN
              ismax=i
              xmax=x(i)
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                          FUNCTION ismin(n,x,inc)
C
C
C***********************************************************************
C
C
C     Function to locate index of minimum element of a real vector.
C
C
C=======================================================================

        INTEGER ismin,n,inc,i
        REAL x(1),xmin

        ismin=1
        xmin=x(1)

        DO 10 i=2,n,inc
           IF(x(i).LT.xmin) THEN
              ismin=i
              xmin=x(i)
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                 FUNCTION isrchigt(n,iarray,inc,itarget)
C
C
C***********************************************************************
C
C
C     Function to return index of first element of iarray greater
C     than itarget, or n+1 if none is found.
C
C
C=======================================================================

        INTEGER isrchigt,n,iarray(1),inc,itarget,i

        isrchigt=n+1

        DO 10 i=1,n,inc
           IF(iarray(i).GT.itarget) THEN
              isrchigt=i
              RETURN
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE wheneq(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).EQ.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenfgt(n,array,inc,target,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of a real
C     vector greater than target.
C
C
C=======================================================================

        INTEGER index(1),n,inc,nval,i
        REAL array(1),target

        nval=0

        DO 10 i=1,n,inc
           IF(array(i).GT.target) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenflt(n,array,inc,target,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of a real
C     vector less than target.
C
C
C=======================================================================

        INTEGER index(1),n,inc,nval,i
        REAL array(1),target

        nval=0

        DO 10 i=1,n,inc
           IF(array(i).LT.target) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenige(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector greater than or equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).GE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenigt(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector greater than itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).GT.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenile(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector less than or equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).LE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenilt(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector less than itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).LT.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenne(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector not equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).NE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END
C***********************************************************************
C
C
                          SUBROUTINE ranset(sd)
C
C
C***********************************************************************
C
C
C     Dummy subroutine to initialize random numbers.
C
C
C=======================================================================

        REAL sd

        RETURN 
        END

C***********************************************************************
C
C
                           FUNCTION ranf(idum)
C
C
C***********************************************************************
C
C
C     Function to return random numbers.
C     copied from the Ran2 routine in Numerical Recipes.
C
C=======================================================================

cc        INTEGER iran
cc        REAL ranf,drand

cc        ranf=DRAND(iran)

cc        RETURN 
cc        END
      implicit real*4 (a-h,o-z), integer*4(i-n)
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RANf=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END


C***********************************************************************
C
C
c                            FUNCTION second()
C
C
C***********************************************************************
C
C
C     Subroutine to return elapsed cpu time.
C
C
C=======================================================================

c        REAL etime,utime,stime,x,second

c        x=etime(utime,stime)

c        second=utime+stime     

c        RETURN 
c        END
C***********************************************************************
C
C
                       SUBROUTINE terror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to terminate the program as the result of a fatal
C     error, close the output files, and dump timing information.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER ierror
        REAL second

C=======================================================================

C   Write error message to the log file and to the terminal.
C   --------------------------------------------------------
        CALL outerror(message)
C            --------
        ierror=-1

        CALL outterm(message,ierror)
C            -------

        CALL outdump(1)
C            -------
        CALL outbods
C            -------

C-----------------------------------------------------------------------
C   Stop timing, output timing data, close files, terminate the
C   simulation.
C-----------------------------------------------------------------------

        cputime1=SECOND()

        CALL outcpu
C            ------
        CALL stopout
C            -------

        STOP
        END
C***********************************************************************
C
C
                      SUBROUTINE outterm(message,n)
C
C
C***********************************************************************
C
C
C     Subroutine to output a message to the terminal and to the
C     terminal emulation file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message
        INTEGER n

C=======================================================================
 
        OPEN(UNIT=utermfil,FILE=termfile,STATUS='OLD')

C   Write the message.
C   ------------------

        IF(n.GE.0) THEN
           IF(.NOT.cosmo) THEN
              WRITE(uterm,*) message,n
              WRITE(utermfil,*) message,n
           ELSE
              WRITE(uterm,*) message,n,' expansion parameter = ',
     &                       expanpar
              WRITE(utermfil,*) message,n,' expansion parameter = ',
     &                          expanpar
           ENDIF
        ELSE
           WRITE(uterm,40)
           WRITE(uterm,50) message 
           WRITE(uterm,40)
           WRITE(utermfil,40)
           WRITE(utermfil,50) message 
           WRITE(utermfil,40)
        ENDIF

 40     FORMAT(/,1x,72('*'))
 50     FORMAT(/,a)

        CLOSE(UNIT=utermfil)

        RETURN
        END
C***********************************************************************
C
C
                  SUBROUTINE outenrgy(am,cmpos,cmvel)
C
C
C***********************************************************************
C
C
C     Subroutine to output diagnostic data to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER k
        REAL am(3),cmpos(ndim),cmvel(ndim),cpunew,cpuold,cpustep,
     &       etotrad,second

        SAVE cpuold

        DATA cpuold/0.0/

C=======================================================================

        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

C-----------------------------------------------------------------------
C   Write mass, energy, angular momentum, center of mass quantities.
C-----------------------------------------------------------------------
 
        etotrad=etot-eradiate-snheat-enerror

        ireclog=ireclog+1
        WRITE(ulog,5,REC=ireclog)

        IF(.NOT.starform) THEN
           ireclog=ireclog+1
           WRITE(ulog,7,REC=ireclog) mtot,mgas
        ELSE
           ireclog=ireclog+1
           WRITE(ulog,10,REC=ireclog) mtot,mgas,mstar
        ENDIF

        ireclog=ireclog+1
        WRITE(ulog,20,REC=ireclog) etotrad,ektot,eptot
        ireclog=ireclog+1
        WRITE(ulog,22,REC=ireclog) ethtot,esofttot

        IF(cosmo.AND.comove) THEN
           ireclog=ireclog+1
           WRITE(ulog,25,REC=ireclog) ecosm
        ENDIF
         
        IF(radiate) THEN
           ireclog=ireclog+1
           WRITE(ulog,30,REC=ireclog) eradiate
        ENDIF

        IF(starform) THEN
           ireclog=ireclog+1
           WRITE(ulog,35,REC=ireclog) snheat
        ENDIF

        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog) am(1),am(2),am(3)
        ireclog=ireclog+1
        WRITE(ulog,50,REC=ireclog) (cmpos(k),k=1,ndim)
        ireclog=ireclog+1
        WRITE(ulog,60,REC=ireclog) (cmvel(k),k=1,ndim)

        IF(usebh) THEN
           ireclog=ireclog+1
           WRITE(ulog,62,REC=ireclog) amvecbh(1),amvecbh(2),amvecbh(3)
           ireclog=ireclog+1
           WRITE(ulog,64,REC=ireclog) (cmposbh(k),k=1,ndim)
           ireclog=ireclog+1
           WRITE(ulog,66,REC=ireclog) (cmvelbh(k),k=1,ndim)
        ENDIF

        cpunew=SECOND()
        cpustep=cpunew-cpuold
        cpuold=cpunew

        ireclog=ireclog+1
        WRITE(ulog,5,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,70,REC=ireclog) cpustep
 
        ireclog=ireclog+1
        WRITE(ulog,5,REC=ireclog)

        CLOSE(UNIT=ulog)

 5      FORMAT(' ')
 7      FORMAT(7x,'mtot,mgas = ',2(1pe12.4))
 10     FORMAT(7x,'mtot,mgas,mstar = ',3(1pe12.4))
 20     FORMAT(7x,'e, ek, ep = ',3(1pe17.9))
 22     FORMAT(7x,'et, es = ',2(1pe17.9))
 25     FORMAT(7x,'ecosm = ',1pe12.4)
 30     FORMAT(7x,'total energy radiated = ',1pe17.9) 
 35     FORMAT(7x,'total supernova heating = ',1pe17.9)
 40     FORMAT(7x,'amx, amy, amz = ',3(1pe17.9))
 50     FORMAT(7x,'cmpos = ',3(1pe17.9))
 60     FORMAT(7x,'cmvel = ',3(1pe17.9))
 62     FORMAT(7x,'amvecbh = ',3(1pe17.9))
 64     FORMAT(7x,'cmposbh = ',3(1pe17.9))
 66     FORMAT(7x,'cmvelbh = ',3(1pe17.9))
 70     FORMAT(10x,'cpu time per step = ',1pe12.4)
 
        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE outerror(message)
C
C
C***********************************************************************
C
C
C     Subroutine to output error messages to the log file.
C
C
C=======================================================================
 
        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*(*) message

C=======================================================================

        OPEN(UNIT=ulog,FILE=logfile,STATUS='OLD',ACCESS='DIRECT',
     &       FORM='FORMATTED',RECL=81)

C   Write the message.
C   ------------------
        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,50,REC=ireclog) message

        ireclog=ireclog+1
        WRITE(ulog,10,REC=ireclog)

        ireclog=ireclog+1
        WRITE(ulog,40,REC=ireclog)
 
        CLOSE(UNIT=ulog)

 10     FORMAT(' ')
 40     FORMAT(1x,72('*'))
 50     FORMAT(a)

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE initent
C
C
C***********************************************************************
C
C
C     Subroutine to initialize the entropy.  Several options are 
C     currently implemented.
C
C     If entinit is .TRUE., then it is assumed that the entropy has
C     not been input and is to be computed using a Jeans mass 
C     criterion.
C
C     If sphinit is .TRUE., then the entropy is smoothed to provide 
C     relatively quiet initial conditions.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),iwsm,iwsm1
        REAL hsminv,wnorm,distnorm,dr2,dr2p,drw,wsm,wmass,dr2i,drw1,
     &       wsm1,wmass1,tcutoff,dx,dy,dz

        EQUIVALENCE (nb(1),bodlist(1))

C=======================================================================

        IF(entinit) THEN

           DO 10 p=1,nsph
              tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(p)*nsmooth/8./ 
     &                (rho(p)/cosmof3))**(2./3.))-((epsvect(p)/  
     &                piinv)**2.)/3.)*(rho(p)/cosmof3) 
              IF(tcutoff.LT.mintemp) tcutoff=mintemp
              entropy(p)=tcutoff/(meanmwt*mhboltz*rho(p)**gamma1)
 10        CONTINUE

        ENDIF

        IF(readinet.AND.(.NOT.entinit)) THEN

           DO 20 p=1,nsph
              entropy(p)=entropy(p)*gamma1/rho(p)**gamma1
 20        CONTINUE

        ELSE

           IF((.NOT.readinas).AND.(.NOT.entinit)) THEN

              DO 30 p=1,nsph
                 entropy(p)=entropy(p)/(meanmwt*mhboltz*
     &                      rho(p)**gamma1)
 30           CONTINUE

           ENDIF

        ENDIF

        IF(sphinit) THEN

C   Zero accumulator thermal energy.
C   --------------------------------
           DO 40 p=1,nsph
              entold(p)=entropy(p)
              entropy(p)=0.
 40        CONTINUE

C   Compute density of particles.
C   -----------------------------

           DO 60 p=1,nsph

              IF(symmetry.EQ.'hk') THEN
                 hsminv=1.0/hsmooth(p)
                 wnorm=piinv*hsminv*hsminv*hsminv
                 distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 50 i=1,nnearlis(p)
                    nb(i)=nearbods(pnear(p)+i-1)
                    dx=pos(p,1)-pos(nb(i),1)
                    dy=pos(p,2)-pos(nb(i),2)
                    dz=pos(p,3)-pos(nb(i),3)
                    IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                    IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                    IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                    IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                    IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                    IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                    dr2=dx**2+dy**2+dz**2
                    dr2p=dr2*distnorm
                    iwsm=dr2p
                    IF(ninterp.LT.iwsm) iwsm=ninterp
                    drw=dr2p-iwsm
                    wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                    wmass=wnorm*wsm
                    dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                    iwsm1=dr2i
                    IF(ninterp.LT.iwsm1) iwsm1=ninterp
                    drw1=dr2i-iwsm1
                    wsm1=(1.-drw1)*wsmooth(iwsm1)+drw1*wsmooth(1+iwsm1)
                    wmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                     hsmooth(nb(i)))*wsm1
                    entropy(p)=entropy(p)+0.5*mass(nb(i))*(wmass+
     &                         wmass1)*entold(nb(i))
                    entropy(nb(i))=entropy(nb(i))+0.5*mass(p)*
     &                             (wmass+wmass1)*entold(p)
 50              CONTINUE

              ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
                 DO 55 i=1,nnearlis(p)
                    nb(i)=nearbods(pnear(p)+i-1)
                    dx=pos(p,1)-pos(nb(i),1)
                    dy=pos(p,2)-pos(nb(i),2)
                    dz=pos(p,3)-pos(nb(i),3)
                    IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                    IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                    IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                    IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                    IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                    IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                    dr2=dx**2+dy**2+dz**2
                    hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                    wnorm=piinv*hsminv*hsminv*hsminv
                    distnorm=hsminv**2*deldr2i
                    dr2p=dr2*distnorm
                    iwsm=dr2p
                    IF(ninterp.LT.iwsm) iwsm=ninterp
                    drw=dr2p-iwsm
                    wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
                    wmass=wnorm*wsm
                    entropy(p)=entropy(p)+mass(nb(i))*wmass*
     &                         entold(nb(i))
                    entropy(nb(i))=entropy(nb(i))+mass(p)*
     &                             wmass*entold(p)
 55              CONTINUE

              ENDIF

C   Include self-interaction term.
C   ------------------------------
              IF(symmetry.EQ.'be') THEN
                 wnorm=piinv/hsmooth(p)**3
              ENDIF

              wmass=wnorm*mass(p)
              entropy(p)=entropy(p)+wmass*entold(p)

 60        CONTINUE

           DO 70 p=1,nsph
              entropy(p)=entropy(p)/rho(p)
 70        CONTINUE

        ENDIF

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE accsphbv
C
C
C***********************************************************************
C
C
C     Subroutine to compute SPH contribution to the acceleration for
C     bulk viscosity.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------
        INTEGER p,i,nb(nbodsmax),nnind(nbodsmax),inear,iwsm,j,
     &          jsubset(nbodsmax),ninear,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,dwsm,dwmass1,aij(nworkvec),
     &       qi,qj,hsmdvvni,vdotdr,hsmdvvp,distnorm,dr2p,dwmass,
     &       dr2i,drw1,dwsm1,aijmass,rhoprhoi,gradprho

        EQUIVALENCE (nb(1),bodlist(1)),(nnind(1),templist(1)),
     &              (jsubset(1),subindex(1)),(aij(1),workvect(1))

C=======================================================================

        DO 10 p=1,nsph
           isubset(p)=0
 10     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 20 p=1,nsphact
           isubset(pactive(p))=1
 20     CONTINUE

C   Compute SPH acceleration, maintaining inversion symmetry.
C   ---------------------------------------------------------
        DO 90 p=1,nsph

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 25 i=1,nnearlis(p)
              jsubset(i)=isubset(nearbods(pnear(p)+i-1))+isubset(p)
 25        CONTINUE

           CALL WHENNE(nnearlis(p),jsubset,1,0,nnind,inear)

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

           IF(symmetry.EQ.'hk') THEN

              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    tempvect(i)=hsmdvvni
                 ELSE
                    tempvect(i)=hsmdvvp
                 ENDIF
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 qi=alpha*csound(p)*hsmdvvp+beta*hsmdvvp**2
                 qj=alpha*csound(nb(i))*hsmdvvni+beta*hsmdvvni**2
                 IF(vdotdr.GT.zero) THEN
                    qi=zero
                    qj=zero
                 ENDIF
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    gradprho=2.*csound(p)*csound(nb(i))/(gamma*
     &                       SQRT(rhoprhoi))
                 ELSE
                    gradprho=csound(p)**2/(gamma*rho(p))+
     &                       csound(nb(i))**2/(gamma*rho(nb(i)))
                 ENDIF
                 aij(i)=(gradprho+(qi/rho(p)+qj/rho(nb(i))))*0.5*
     &                  (dwmass+dwmass1)
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 40 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    tempvect(i)=hsmdvvni
                 ELSE
                    tempvect(i)=hsmdvvp
                 ENDIF
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 qi=alpha*csound(p)*hsmdvvp+beta*hsmdvvp**2
                 qj=alpha*csound(nb(i))*hsmdvvni+beta*hsmdvvni**2
                 IF(vdotdr.GT.zero) THEN
                    qi=zero
                    qj=zero
                 ENDIF
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    gradprho=2.*csound(p)*csound(nb(i))/(gamma*
     &                       SQRT(rhoprhoi))
                 ELSE
                    gradprho=csound(p)**2/(gamma*rho(p))+
     &                       csound(nb(i))**2/(gamma*rho(nb(i)))
                 ENDIF
                 aij(i)=(gradprho+(qi/rho(p)+qj/rho(nb(i))))*dwmass
 40           CONTINUE

           ENDIF

           IF(isubset(p).EQ.1) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 50 i=1,inear
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 aijmass=aij(i)*mass(nb(i))*cosmofac
                 acc(p,1)=acc(p,1)-aijmass*dx
                 acc(p,2)=acc(p,2)-aijmass*dy
                 acc(p,3)=acc(p,3)-aijmass*dz
 50           CONTINUE
           ENDIF

           CALL WHENEQ(inear,jsubset,1,1,nnind,ninear)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,ninear
              j=nnind(i)
              dx=pos(p,1)-pos(nb(j),1)
              dy=pos(p,2)-pos(nb(j),2)
              dz=pos(p,3)-pos(nb(j),3)
              IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
              IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
              IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
              IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
              IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
              IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
              aijmass=aij(j)*mass(p)*cosmofac
              acc(nb(j),1)=acc(nb(j),1)+aijmass*dx
              acc(nb(j),2)=acc(nb(j),2)+aijmass*dy
              acc(nb(j),3)=acc(nb(j),3)+aijmass*dz
 60        CONTINUE

 90     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE ethdotbv(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the time derivative of the thermal
C     energy for bulk viscosity.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,nb(nbodsmax),iwsm,imumax,ismax,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,dwsm,dwsm1,qi,qj,vdotdr,
     &       tmuij,hsmdvvni,dr2p,hsmdvvp,distnorm,dwmass,dr2i,drw1,
     &       dwmass1,eijp,eijnbi,rhoprhoi

        EQUIVALENCE (nb(1),bodlist(1))

C=======================================================================

C   Zero out derivative of thermal energy.
C   --------------------------------------
        DO 10 p=1,nsph
           dethdt(p)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute derivative of thermal energy
C-----------------------------------------------------------------------

        DO 40 p=1,nsph

           IF(symmetry.EQ.'hk') THEN
              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    tempvect(i)=hsmdvvni
                 ELSE
                    tempvect(i)=hsmdvvp
                 ENDIF
                 qi=alpha*csound(p)*hsmdvvp+beta*hsmdvvp**2
                 qj=alpha*csound(nb(i))*hsmdvvni+beta*
     &              hsmdvvni**2
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 IF(vdotdr.GT.zero) THEN
                    qi=zero
                    qj=zero
                 ENDIF
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    eijp=vdotdr*(csound(nb(i))*csound(p)/
     &                   (SQRT(rhoprhoi)*gamma)+0.5*(qi/rho(p)+
     &                   qj/rho(nb(i))))
                    eijnbi=eijp
                 ELSE
                    eijp=vdotdr*(csound(p)**2/(gamma*rho(p))+
     &                   0.5*(qi/rho(p)+qj/rho(nb(i))))
                    eijnbi=vdotdr*(csound(nb(i))**2/(gamma*rho(nb(i)))+
     &                   0.5*(qi/rho(p)+qj/rho(nb(i))))
                 ENDIF
                 dethdt(p)=dethdt(p)+eijp*mass(nb(i))*0.5*(dwmass+
     &                     dwmass1)
                 dethdt(nb(i))=dethdt(nb(i))+eijnbi*mass(p)*0.5*(dwmass+
     &                         dwmass1)
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    tempvect(i)=hsmdvvni
                 ELSE
                    tempvect(i)=hsmdvvp
                 ENDIF
                 qi=alpha*csound(p)*hsmdvvp+beta*hsmdvvp**2
                 qj=alpha*csound(nb(i))*hsmdvvni+beta*
     &              hsmdvvni**2
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 IF(vdotdr.GT.zero) THEN
                    qi=zero
                    qj=zero
                 ENDIF
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    eijp=vdotdr*(csound(nb(i))*csound(p)/
     &                   (SQRT(rhoprhoi)*gamma)+0.5*(qi/rho(p)+
     &                   qj/rho(nb(i))))
                    eijnbi=eijp
                 ELSE
                    eijp=vdotdr*(csound(p)**2/(gamma*rho(p))+
     &                   0.5*(qi/rho(p)+qj/rho(nb(i))))
                    eijnbi=vdotdr*(csound(nb(i))**2/(gamma*rho(nb(i)))+
     &                   0.5*(qi/rho(p)+qj/rho(nb(i))))
                 ENDIF
                 dethdt(p)=dethdt(p)+eijp*mass(nb(i))*dwmass
                 dethdt(nb(i))=dethdt(nb(i))+eijnbi*mass(p)*dwmass
 35           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0.AND.pc.EQ.'correct') THEN
              imumax=ISMAX(nnearlis(p),tempvect,1)
              tmuij=tempvect(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 60 i=1,nnearlis(p)
                 tmuij=tempvect(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 60           CONTINUE

           ENDIF

 40     CONTINUE

C##### LH 050990
                 
        DO 70 p=1,nsph
           dethdt(p)=dethdt(p)/cosmofac
 70     CONTINUE
				                         
C##### LH 050990

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE accsphcv
C
C
C***********************************************************************
C
C
C     Subroutine to compute SPH contribution to the acceleration.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        INTEGER p,i,nb(nbodsmax),nnind(nbodsmax),inear,iwsm,
     &          j,jsubset(nbodsmax),ninear,iwsm1
        REAL dx,dy,dz,hsmavg,dr2,hsminv,drw,dwnorm,dwsm1,muij(nbodsmax),
     &       dwsm,vdotdr,aij(nworkvec),piij,distnorm,dr2p,dwmass,
     &       dr2i,drw1,dwmass1,aijmass,formp,formnbi,abshdivv,rhoprhoi,
     &       gradprho

        EQUIVALENCE (nb(1),bodlist(1)),(nnind(1),templist(1)),
     &              (jsubset(1),subindex(1)),(muij(1),tempvect(1)),
     &              (aij(1),workvect(1))

C=======================================================================

        DO 10 p=1,nsph
           isubset(p)=0
 10     CONTINUE

CVD$ NODEPCHK
CDIR$ IVDEP
        DO 20 p=1,nsphact
           isubset(pactive(p))=1
 20     CONTINUE

C   Compute SPH acceleration, maintaining inversion symmetry.
C   ---------------------------------------------------------
        DO 90 p=1,nsph

           abshdivv=ABS(hsmdivv(p))
           formp=abshdivv/(abshdivv+hsmcurlv(p)+0.0001*csound(p))

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 25 i=1,nnearlis(p)
              jsubset(i)=isubset(nearbods(pnear(p)+i-1))+isubset(p)
 25        CONTINUE

           CALL WHENNE(nnearlis(p),jsubset,1,0,nnind,inear)

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

           IF(symmetry.EQ.'hk') THEN

              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/
     &                   (dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    gradprho=2.*csound(p)*csound(nb(i))/(gamma*
     &                       SQRT(rhoprhoi))
                 ELSE
                    gradprho=csound(p)**2/(gamma*rho(p))+
     &                       csound(nb(i))**2/(gamma*rho(nb(i)))
                 ENDIF
                 aij(i)=(gradprho+piij)*0.5*(dwmass+dwmass1)
                 muij(i)=ABS(muij(i))
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 40 i=1,inear
                 nb(i)=nearbods(pnear(p)+nnind(i)-1)
                 jsubset(i)=isubset(nb(i))
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/
     &                   (dr2+epssph*hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    gradprho=2.*csound(p)*csound(nb(i))/(gamma*
     &                       SQRT(rhoprhoi))
                 ELSE
                    gradprho=csound(p)**2/(gamma*rho(p))+
     &                       csound(nb(i))**2/(gamma*rho(nb(i)))
                 ENDIF
                 aij(i)=(gradprho+piij)*dwmass
                 muij(i)=ABS(muij(i))
 40           CONTINUE

           ENDIF

           IF(isubset(p).EQ.1) THEN

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 50 i=1,inear
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 aijmass=aij(i)*mass(nb(i))*cosmofac
                 acc(p,1)=acc(p,1)-aijmass*dx
                 acc(p,2)=acc(p,2)-aijmass*dy
                 acc(p,3)=acc(p,3)-aijmass*dz
 50           CONTINUE
           ENDIF

           CALL WHENEQ(inear,jsubset,1,1,nnind,ninear)

CVD$ NODEPCHK
CDIR$ IVDEP
           DO 60 i=1,ninear
              j=nnind(i)
              dx=pos(p,1)-pos(nb(j),1)
              dy=pos(p,2)-pos(nb(j),2)
              dz=pos(p,3)-pos(nb(j),3)
              IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
              IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
              IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
              IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
              IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
              IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
              aijmass=aij(j)*mass(p)*cosmofac
              acc(nb(j),1)=acc(nb(j),1)+aijmass*dx
              acc(nb(j),2)=acc(nb(j),2)+aijmass*dy
              acc(nb(j),3)=acc(nb(j),3)+aijmass*dz
 60        CONTINUE

 90     CONTINUE

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE ethdotcv(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the time derivative of the thermal
C     energy.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,nb(nbodsmax),iwsm,imumax,ismax,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,muij(nbodsmax),dwsm,vdotdr,
     &       piij,tmuij,dr2p,distnorm,dwmass,dr2i,drw1,dwsm1,dwmass1,
     &       hsmavg,formp,formnbi,abshdivv,eijp,eijnbi,rhoprhoi

        EQUIVALENCE (nb(1),bodlist(1)),(muij(1),tempvect(1))

C=======================================================================

C   Zero out derivative of thermal energy.
C   --------------------------------------
        DO 10 p=1,nsph
           dethdt(p)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute derivative of thermal energy.
C-----------------------------------------------------------------------

        DO 40 p=1,nsph

           abshdivv=ABS(hsmdivv(p))
           formp=abshdivv/(abshdivv+hsmcurlv(p)+0.0001*csound(p))

           IF(symmetry.EQ.'hk') THEN
              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/(dr2+epssph*
     &                   hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    eijp=vdotdr*(csound(nb(i))*csound(p)/
     &                   (SQRT(rhoprhoi)*gamma)+0.5*piij)
                    eijnbi=eijp
                 ELSE
                    eijp=vdotdr*(csound(p)**2/(gamma*rho(p))+0.5*piij)
                    eijnbi=vdotdr*(csound(nb(i))**2/(gamma*rho(nb(i)))+
     &                   0.5*piij)
                 ENDIF
                 dethdt(p)=dethdt(p)+eijp*mass(nb(i))*0.5*(dwmass+
     &                     dwmass1)
                 dethdt(nb(i))=dethdt(nb(i))+eijnbi*mass(p)*0.5*(dwmass+
     &                         dwmass1)
                 muij(i)=ABS(muij(i))
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/(dr2+epssph*
     &                   hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 IF(geogradp) THEN
                    rhoprhoi=rho(p)*rho(nb(i))
                    eijp=vdotdr*(csound(nb(i))*csound(p)/
     &                   (SQRT(rhoprhoi)*gamma)+0.5*piij)
                    eijnbi=eijp
                 ELSE
                    eijp=vdotdr*(csound(p)**2/(gamma*rho(p))+0.5*piij)
                    eijnbi=vdotdr*(csound(nb(i))**2/(gamma*rho(nb(i)))+
     &                   0.5*piij)
                 ENDIF
                 dethdt(p)=dethdt(p)+eijp*mass(nb(i))*dwmass
                 dethdt(nb(i))=dethdt(nb(i))+eijnbi*mass(p)*dwmass
                 muij(i)=ABS(muij(i))
 35           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0.AND.pc.EQ.'correct') THEN
              imumax=ISMAX(nnearlis(p),muij,1)
              tmuij=muij(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 60 i=1,nnearlis(p)
                 tmuij=muij(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 60           CONTINUE

           ENDIF

 40     CONTINUE

C##### LH 050990
                 
        DO 70 p=1,nsph
           dethdt(p)=dethdt(p)/cosmofac
 70     CONTINUE
				                         
C##### LH 050990

        RETURN
        END


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C   Description:
C   -----------
C 
C   This routine finds the cubic spline coefficients.  The boundary condi-
C   tions may be one of the following three:
C        1) "natural," that is, zero second derivatives
C        2) first derivatives specified
C        3) third derivatives computed from supplied data
C 
C 
C   Call sequence:
C   -------------
C   call spline(x,y,n,yp1,ypn,y2)
C 
C   integer n
C   real x(n),y(n),yp1,ypn,y2(n)
C 
C   Parameters:
C   ----------
C 
C   n        number of supplied grid points
C   x        abcissa array
C   y        ordinate array
C   yp1      boundary condition at j=1
C   ypn      boundary condition at j=n
C   y2       array to contain spline coefficients
C 
C   Returns:
C   -------
C 
C   None, spline coefficients returned by pointer
C 
C   Notes:
C   -----
C 
C   If     yp1,yp2 >  1.0e30  boundary conditions (1) natural splines are used
C   If     yp1,yp2 < -1.0e30  boundary conditions (3) approx. 3rd derivs used
C   Otherwise                 boundary conditions (2) explicit 2nd derivs used
C 
C   By:
C   --
C   Adopted from Numerical Recipes, Press et al.
C   Third deriv. boundary condition ---  MDW 11/13/88
C 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C          Routines to compute spline coefficients interpolate and
C          integrate function
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE                SPLINE(X,Y,N,YP1,YPN,Y2)
C
C          Variables:
C          ---------
C          X          Input mesh array
C          Y          Input function values
C          N          Dimension of arrays
C          YP1        First derivative at X(1)
C          YPN        First derivative at X(N)
C          Y2         Second derivatives (spline coefficients)
C
      IMPLICIT                  LOGICAL (A - Z)
      INTEGER                   N,NMAX,I,K
      PARAMETER                 (NMAX=1000)
      REAL                      X(N),Y(N),Y2(N),U(NMAX),YP1,YPN,SIG,P,
     &                          QN,UN,D1,D2
C
      IF (YP1.LT.-0.99D20) THEN
         Y2(1)=1.0
         D2 = ((Y(4)-Y(3))/(X(4)-X(3)) - (Y(3)-Y(2))/
     &        (X(3)-X(2)))/(X(4)-X(2))
         D1 = ((Y(3)-Y(2))/(X(3)-X(2)) - (Y(2)-Y(1))/
     &        (X(2)-X(1)))/(X(3)-X(1))
         U(1) = -6.0*(D2-D1)*(X(2)-X(1))/(X(4)-X(1))
      ELSE IF (YP1.GT..99D20) THEN
           Y2(1)=0.0D0
           U(1)=0.0D0
      ELSE
           Y2(1)=-0.5D0
           U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      END IF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.D0
        Y2(I)=(SIG-1.D0)/P
        U(I)=(6.D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99D20) THEN
           QN=0.0D0
           UN=0.0D0
      ELSE
           QN=0.5D0
           UN=(3.D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      END IF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.D0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
C***********************************************************************
C
C
                          SUBROUTINE entdotbv(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the time derivative of the entropic
C     function a(s) using the bulk viscosity.
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,nb(nbodsmax),iwsm,imumax,ismax,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,dwsm,dwsm1,qi,qj,vdotdr,
     &       tmuij,hsmdvvni,dr2p,hsmdvvp,distnorm,dwmass,dr2i,drw1,
     &       dwmass1,eij

        EQUIVALENCE (nb(1),bodlist(1))

C=======================================================================

C   Zero out derivative of thermal energy.
C   --------------------------------------
        DO 10 p=1,nsph
           dentdt(p)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute derivative of thermal energy
C-----------------------------------------------------------------------

        DO 40 p=1,nsph

           IF(symmetry.EQ.'hk') THEN
              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    tempvect(i)=hsmdvvni
                 ELSE
                    tempvect(i)=hsmdvvp
                 ENDIF
                 qi=alpha*csound(p)*hsmdvvp+beta*hsmdvvp**2
                 qj=alpha*csound(nb(i))*hsmdvvni+beta*
     &              hsmdvvni**2
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 IF(vdotdr.GT.zero) THEN
                    qi=zero
                    qj=zero
                 ENDIF
                 eij=rho(p)*rho(nb(i))
                 eij=vdotdr*0.5*(qi/rho(p)+qj/rho(nb(i)))
                 dentdt(p)=dentdt(p)+eij*mass(nb(i))*0.5*(dwmass+
     &                     dwmass1)
                 dentdt(nb(i))=dentdt(nb(i))+eij*mass(p)*0.5*(dwmass+
     &                         dwmass1)
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 hsmdvvni=hsmdivv(nb(i))
                 IF(hsmdvvni.GT.zero) THEN
                    hsmdvvni=zero
                 ELSE
                    hsmdvvni= - hsmdvvni
                 ENDIF
                 IF(hsmdivv(p).GT.zero) THEN
                    hsmdvvp=zero
                 ELSE
                    hsmdvvp= - hsmdivv(p)
                 ENDIF
                 IF(hsmdvvni.GT.hsmdvvp) THEN
                    tempvect(i)=hsmdvvni
                 ELSE
                    tempvect(i)=hsmdvvp
                 ENDIF
                 qi=alpha*csound(p)*hsmdvvp+beta*hsmdvvp**2
                 qj=alpha*csound(nb(i))*hsmdvvni+beta*
     &              hsmdvvni**2
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 IF(vdotdr.GT.zero) THEN
                    qi=zero
                    qj=zero
                 ENDIF
                 eij=rho(p)*rho(nb(i))
                 eij=vdotdr*0.5*(qi/rho(p)+qj/rho(nb(i)))
                 dentdt(p)=dentdt(p)+eij*mass(nb(i))*dwmass
                 dentdt(nb(i))=dentdt(nb(i))+eij*mass(p)*dwmass
 35           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0.AND.pc.EQ.'correct') THEN
              imumax=ISMAX(nnearlis(p),tempvect,1)
              tmuij=tempvect(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 60 i=1,nnearlis(p)
                 tmuij=tempvect(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 60           CONTINUE

           ENDIF

 40     CONTINUE

C##### LH 050990

        DO 70 p=1,nsph
           dentdt(p)=dentdt(p)*gamma1/(cosmofac*rho(p)**gamma1)
 70     CONTINUE

C##### LH 050990

        RETURN
        END
C***********************************************************************
C
C
                          SUBROUTINE entdotcv(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to compute the time derivative of the entropic
C     function a(s).
C
C
C=======================================================================

        INCLUDE 'treedefs.h'

C   Declaration of local variables.
C   -------------------------------

        CHARACTER*7 pc
        INTEGER p,i,nb(nbodsmax),iwsm,imumax,ismax,iwsm1
        REAL dx,dy,dz,dr2,hsminv,drw,dwnorm,muij(nbodsmax),dwsm,vdotdr,
     &       piij,tmuij,dr2p,distnorm,dwmass,dr2i,drw1,dwsm1,dwmass1,
     &       hsmavg,eij,formp,formnbi,abshdivv

        EQUIVALENCE (nb(1),bodlist(1)),(muij(1),tempvect(1))

C=======================================================================

C   Zero out derivative of thermal energy.
C   --------------------------------------
        DO 10 p=1,nsph
           dentdt(p)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute derivative of thermal energy.
C-----------------------------------------------------------------------

        DO 40 p=1,nsph

           abshdivv=ABS(hsmdivv(p))
           formp=abshdivv/(abshdivv+hsmcurlv(p)+0.0001*csound(p))

           IF(symmetry.EQ.'hk') THEN

              hsminv=1.0/hsmooth(p)
              dwnorm=piinv*hsminv*hsminv*hsminv*hsminv*hsminv
              distnorm=hsminv*hsminv*deldr2i

C   Compute interaction between body p and neighbors.
C   -------------------------------------------------

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 30 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 dr2i=dr2*deldr2i/(hsmooth(nb(i))*hsmooth(nb(i)))
                 iwsm1=dr2i
                 IF(ninterp.LT.iwsm1) iwsm1=ninterp
                 drw1=dr2i-iwsm1
                 dwsm1=(1.-drw1)*dwsmooth(iwsm1)+drw1*dwsmooth(1+iwsm1)
                 dwmass1=piinv/(hsmooth(nb(i))*hsmooth(nb(i))*
     &                   hsmooth(nb(i))*hsmooth(nb(i))*hsmooth(nb(i)))*
     &                   dwsm1
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/(dr2+epssph*
     &                   hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 eij=rho(p)*rho(nb(i))
                 eij=vdotdr*0.5*piij
                 dentdt(p)=dentdt(p)+eij*mass(nb(i))*0.5*(dwmass+
     &                     dwmass1)
                 dentdt(nb(i))=dentdt(nb(i))+eij*mass(p)*0.5*(dwmass+
     &                         dwmass1)
                 muij(i)=ABS(muij(i))
 30           CONTINUE

           ELSE

CVD$ NODEPCHK
CDIR$ IVDEP
              DO 35 i=1,nnearlis(p)
                 nb(i)=nearbods(pnear(p)+i-1)
                 dx=pos(p,1)-pos(nb(i),1)
                 dy=pos(p,2)-pos(nb(i),2)
                 dz=pos(p,3)-pos(nb(i),3)
                 IF(dx.GE.hboxsize.AND.bwrap) dx=dx-pboxsize
                 IF(dx.LT.-hboxsize.AND.bwrap) dx=dx+pboxsize
                 IF(dy.GE.hboxsize.AND.bwrap) dy=dy-pboxsize
                 IF(dy.LT.-hboxsize.AND.bwrap) dy=dy+pboxsize
                 IF(dz.GE.hboxsize.AND.bwrap) dz=dz-pboxsize
                 IF(dz.LT.-hboxsize.AND.bwrap) dz=dz+pboxsize
                 dr2=dx*dx+dy*dy+dz*dz
                 hsminv=two/(hsmooth(p)+hsmooth(nb(i)))
                 dwnorm=piinv*hsminv**2*hsminv**2*hsminv
                 distnorm=hsminv**2*deldr2i
                 dr2p=dr2*distnorm
                 iwsm=dr2p
                 IF(ninterp.LT.iwsm) iwsm=ninterp
                 drw=dr2p-iwsm
                 dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
                 dwmass=dwnorm*dwsm
                 vdotdr=(veltpos(p,1)-veltpos(nb(i),1))*dx+
     &                  (veltpos(p,2)-veltpos(nb(i),2))*dy+
     &                  (veltpos(p,3)-veltpos(nb(i),3))*dz+
     &                  (dx*dx+dy*dy+dz*dz)*cosmohub
                 abshdivv=hsmdivv(nb(i))
                 abshdivv=ABS(abshdivv)
                 formnbi=abshdivv/(abshdivv+hsmcurlv(nb(i))+0.0001*
     &                   csound(nb(i)))
                 hsmavg=0.5*(hsmooth(p)+hsmooth(nb(i)))
                 muij(i)=vdotdr*hsmavg*0.5*(formp+formnbi)/(dr2+epssph*
     &                   hsmavg**2)
                 IF(vdotdr.GT.zero) muij(i)=zero
                 piij=(-alpha*muij(i)*(0.5*(csound(p)+csound(nb(i))))+
     &                beta*muij(i)**2)/(0.5*(rho(p)+rho(nb(i))))
                 eij=rho(p)*rho(nb(i))
                 eij=vdotdr*0.5*piij
                 dentdt(p)=dentdt(p)+eij*mass(nb(i))*dwmass
                 dentdt(nb(i))=dentdt(nb(i))+eij*mass(p)*dwmass
                 muij(i)=ABS(muij(i))
 35           CONTINUE

           ENDIF

           IF(nnearlis(p).GT.0.AND.pc.EQ.'correct') THEN
              imumax=ISMAX(nnearlis(p),muij,1)
              tmuij=muij(imumax)
              IF(mumaxdvh(p).LT.tmuij) mumaxdvh(p)=tmuij
CVD$ NODEPCHK
CDIR$ IVDEP
              DO 60 i=1,nnearlis(p)
                 tmuij=muij(i)
                 IF(mumaxdvh(nb(i)).LT.tmuij) mumaxdvh(nb(i))=tmuij
 60           CONTINUE

           ENDIF

 40     CONTINUE

C##### LH 050990

        DO 70 p=1,nsph
           dentdt(p)=dentdt(p)*gamma1/(cosmofac*rho(p)**gamma1)
 70     CONTINUE

C##### LH 050990

        RETURN
        END
C***********************************************************************
C
C
                       SUBROUTINE stepent(pc)
C
C
C***********************************************************************
C
C
C     Subroutine to advance the entropy implicitly by a timestep dt, 
C     using the Newton-Raphson method.  The local variable p is a 
C     pointer to the bodies.
C
C     In this version of stepent, and in the subroutines ethrad and 
C     entdiff (the latter called by bisectas), radiative cooling is 
C     cut off below a critical temperature, tcutoff, determined by a 
C     Jeans mass criterion.  Specifically, if the Jeans mass is 
C     written in the form:
C
C        M-J = mjconst * csound**4 / SQRT(P * G**3)
C
C     then for a given M-J the cutoff temperature is given by:
C
C        tcutoff = ctcutoff * meanmwt * mhboltz * G * 
C                  (M-J)**(2/3) * rho **(1/3)
C
C     where meanmwt and mhboltz are as defined in the main routine,
C     G = 1 in the adopted system of units.  For the classical 
C     definition of the Jeans mass, mjconst = pi **(3/2) while for 
C     the Bonnor-Ebert Jeans mass, mjconst = 1.18.  The constant 
C     ctcutoff, input via the parameter file, includes not only 
C     (mjconst)**(-2/3), but a factor determining the number of 
C     particles to include in the Jeans mass.  For the two choices 
C     above, (mjconst)**(-2/3) = 0.31831 and 0.89553, respectively.
C     Typically, the factor from the Jeans mass will be ~ 1 --> 
C     (nsmooth/8) ** (2/3), implying a total variation in ctcutoff 
C     of ~ 10.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        CHARACTER*7 pc
        INTEGER p,nit,i,ib,nitmax,nsubset,nblist,ifirst
        REAL dt,smallnum,temp,templog,tenf1,tenf2,heatcool,dheatc,
     &       tempmin,tempmax,entlast,ethmin,ethmax,enttemp,tempdiff,
     &       cutoff,tendiff,tenf3,tenf3exp,compcool,aconst,tcutoff,
     &       tcutoffc,compdiff,comdiff,cutoffc,rhot
     
        PARAMETER(nitmax=15)
     
        SAVE smallnum,tempmin,tempmax,ifirst,ethmin,ethmax
     
        DATA smallnum/1.e-10/,tempmin/10.0/,tempmax/1.e10/,ifirst/1/
     
C=======================================================================
     
        dt=tnow-teth
     
        IF(ifirst.EQ.1) THEN
           ethmin=tempmin/(gamma1*meanmwt*mhboltz)
           ethmax=tempmax/(gamma1*meanmwt*mhboltz)
           ifirst=0
        ENDIF
     
        IF(pc.EQ.'predict') THEN
           DO 5 p=1,nsph
              dentold(p)=dentdt(p)
 5         CONTINUE
        ENDIF

        DO 10 p=1,nsph
           entold(p)=entropy(p)
           bodlist(p)=p
 10     CONTINUE

        nit=0
        nblist=nsph
     
 20     CONTINUE
     
        IF(.NOT.radiate) THEN
     
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 30 i=1,nblist
              ib=bodlist(i)
              entropy(ib)=entold(ib)+0.5*dt*(dentdt(ib)+dentold(ib))
              enttemp=entropy(ib)
              IF(enttemp.LE.zero) THEN
                 tempvect(i)=two*smallnum
              ELSE
                 tempvect(i)=zero
              ENDIF
 30        CONTINUE
     
        ELSE
     
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 40 i=1,nblist
              ib=bodlist(i)
              temp=meanmwt*mhboltz*entropy(ib)*rho(ib)**gamma1
              tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(ib)*nsmooth/8./
     &                (rho(ib)/cosmof3))**(2./3.))-((epsvect(ib)/
     &                piinv)**2.)/3.)*(rho(ib)/cosmof3)
              IF(tcutoff.LT.mintemp) tcutoff=mintemp
              templog=LOG10(temp)
              tempdiff=10.0*(temp-tcutoff)/(tcutoff+tiny)
              IF(tempdiff.GT.thirty5) tempdiff=thirty5
              IF(tempdiff.LT.-thirty5) tempdiff = -thirty5
              tendiff=10.**tempdiff
              cutoff=0.5*(tendiff-1.0)/(tendiff+1.0)+0.5

              tcutoffc=tcutoff
              IF(tcutoff.LT.10000.) tcutoffc=10000.
              compdiff=10.0*(temp-tcutoffc)/(tcutoffc+tiny)
              IF(compdiff.GT.thirty5) compdiff=thirty5
              IF(compdiff.LT.-thirty5) compdiff=-thirty5
              comdiff=10.**compdiff
              cutoffc=0.5*(comdiff-1.0)/(comdiff+1.0)+0.5

              tenf1=10.**(-4.43-0.273*(4.-templog)**2)
              tenf2=10.**(-0.10-1.880*(5.23-templog)**4)
              tenf3exp=6.2-templog
              IF(templog-6.2.GT.zero) tenf3exp=zero
              tenf3=10.**(-2.*tenf3exp**4-1.7)
              compcool=comptmh*((1.+redshift)**4)

C##### LH 050990

              heatcool=aheatmh*fhydrogn+bheatmh2*rho(ib)*fhydrog2/
     &                 cosmof3-ccoolmh2*rho(ib)*(tenf1+tenf2+tenf3)*
     &                 cutoff*fhydrog2/cosmof3-compcool*temp*cutoffc

C##### LH 050990

              aconst=slowcool*rho(ib)**gamma1*entropy(ib)/(gamma1*dt)+
     &               rho(ib)**gamma1*dentdt(ib)/gamma1
              IF(aconst.LE.zero) aconst=tiny
              aconst=heatcool/aconst
              heatcool=heatcool/SQRT(1.+aconst*aconst)
              tempvect(i)=entropy(ib)-entold(ib)-0.5*dt*(dentdt(ib)+
     &                    gamma1*heatcool/rho(ib)**gamma1+dentold(ib)+
     &                    gamma1*derad(ib)/rho(ib)**gamma1)

C##### LH 050990

              dheatc= - ccoolmh2*rho(ib)/cosmof3*(((0.546*(4.-
     &                templog)*tenf1+
     &                7.52*((5.23-templog)**3)*tenf2+8.*tenf3exp**3*
     &                tenf3)*cutoff*gamma1/(rho(ib)**gamma1*
     &                entropy(ib)))+((tenf1+tenf2+tenf3)*(tendiff/
     &                ((tendiff+1.)**2))*10.0*2.30258509*meanmwt*
     &                mhboltz*gamma1/tcutoff))*fhydrog2-compcool*
     &                meanmwt*mhboltz*gamma1*
     &                (cutoffc+compcool*temp*(comdiff/((comdiff+1.)**2))
     &                *10.0*2.30258509/tcutoffc)

C##### LH 050990

              dheatc=(dheatc*(1.-aconst*aconst/(1.+aconst*aconst))+
     &                aconst*aconst*aconst/(1.+aconst*aconst)*slowcool
     &                /dt)/SQRT(1.+aconst*aconst)
              tempvect(i)=tempvect(i)/(1.-0.5*dt*dheatc)
 40        CONTINUE
     
CVD$ NODEPCHK
CDIR$ IVDEP
           DO 50 i=1,nblist
              ib=bodlist(i)
              entlast=entropy(ib)
              entropy(ib)=entropy(ib)-tempvect(i)
              enttemp=entropy(ib)
              IF(ethmax.LT.rho(ib)**gamma1*enttemp/gamma1) 
     &           enttemp=gamma1*ethmax/rho(ib)**gamma1
              IF(rho(ib)**gamma1*enttemp/gamma1.LT.ethmin) 
     &           enttemp=gamma1*ethmin/rho(ib)**gamma1
              entropy(ib)=enttemp
              tempvect(i)=(entlast-entropy(ib))/entlast
              tempvect(i)=ABS(tempvect(i))
              IF(rho(ib)**gamma1*enttemp/gamma1.LT.1.05*ethmin.OR.
     &           enttemp.LT.entold(ib)/10.) tempvect(i)=two*smallnum
              IF(rho(ib)**gamma1*enttemp/gamma1.GT.0.95*ethmax.OR.
     &           enttemp.GT.10.*entold(ib)) tempvect(i)=two*smallnum
 50        CONTINUE
     
        ENDIF
     
        nit = nit + 1
     
        CALL WHENFGT(nblist,tempvect,1,smallnum,isubset,nsubset)
     
        DO 60 i=1,nsubset
           templist(i)=bodlist(isubset(i))
 60     CONTINUE
     
        DO 70 i=1,nsubset
           bodlist(i)=templist(i)
 70     CONTINUE
     
        nblist=nsubset
     
        IF(nit.LE.nitmax) THEN
     
           IF(nblist.GT.0) GO TO 20
     
        ELSE

           CALL bisectas(nblist,tempmin,tempmax,smallnum,dt)
C               --------     
     
        ENDIF
     
        DO 75 p=1,nsph
           IF(radiate) CALL entcheck(p,dt,smallnum)
C                           --------
 75     CONTINUE

C   Recompute sound speed, pressure / rho**2.
C   -----------------------------------------
        DO 80 p=1,nsph
           rhot=rho(p)*cosmof3
           csound(p)=SQRT(gamma*rhot**gamma1*entropy(p))
 80     CONTINUE
     
        IF(pc.EQ.'correct') THEN
           teth=teth+dt
        ENDIF

        IF(pc.EQ.'predict') THEN
           DO 90 p=1,nsph
              entropy(p)=entold(p)
 90        CONTINUE
        ENDIF
     
        RETURN
        END
C***********************************************************************
C
C
                     SUBROUTINE entcheck(ib,dt,smallnum)
C
C
C***********************************************************************
C
C
C     Subroutine to check whether the new value for the entropy is
C     actually a root.
C
C
C=======================================================================
     
        INCLUDE 'treedefs.h'
     
C   Declaration of local variables.
C   -------------------------------
     
        INTEGER ib
        REAL dt,tenf1,tenf2,heatcool,tempdiff,cutoff,tendiff,tenf3,
     &       tenf3exp,temp,logtemp,flogt,compcool,aconst,tcutoff,
     &       smallnum,tcutoffc,compdiff,comdiff,cutoffc

C=======================================================================

        temp=entropy(ib)*meanmwt*mhboltz*rho(ib)**gamma1
        logtemp=LOG10(temp)

        tcutoff=meanmwt*mhboltz*ctcutoff*(((mass(ib)*nsmooth/8./
     &          (rho(ib)/cosmof3))**(2./3.))-((epsvect(ib)/
     &          piinv)**2.)/3.)*(rho(ib)/cosmof3)
        IF(tcutoff.LT.mintemp) tcutoff=mintemp
        tempdiff=10.0*(temp-tcutoff)/(tcutoff+tiny)
        IF(tempdiff.GT.thirty5) tempdiff=thirty5
        IF(tempdiff.LT.-thirty5) tempdiff = -thirty5
        tendiff=10.**tempdiff
        cutoff=0.5*(tendiff-1.0)/(tendiff+1.0)+0.5

        tcutoffc=tcutoff
        IF(tcutoff.LT.10000.) tcutoffc=10000.
        compdiff=10.0*(temp-tcutoffc)/(tcutoffc+tiny)
        IF(compdiff.GT.thirty5) compdiff=thirty5
        IF(compdiff.LT.-thirty5) compdiff=-thirty5
        comdiff=10.**compdiff
        cutoffc=0.5*(comdiff-1.0)/(comdiff+1.0)+0.5

        tenf1=10.**(-4.43-0.273*(4.-logtemp)**2)
        tenf2=10.**(-0.10-1.880*(5.23-logtemp)**4)
        tenf3exp=6.2-logtemp
        IF(logtemp-6.2.GT.zero) tenf3exp=zero
        tenf3=10.**(-2.*tenf3exp**4-1.7)
        compcool=comptmh*temp*((1.+redshift)**4)

C##### LH 050990

        heatcool=aheatmh*fhydrogn+bheatmh2*rho(ib)*fhydrog2/
     &           cosmof3-ccoolmh2*rho(ib)*(tenf1+tenf2+tenf3)*
     &           cutoff*fhydrog2/cosmof3-compcool*cutoffc

C##### LH 050990

        aconst=slowcool*temp/(meanmwt*mhboltz*gamma1*dt)+
     &         rho(ib)**gamma1*dentdt(ib)/gamma1
        IF(aconst.LE.zero) aconst=tiny
        aconst=heatcool/aconst
        heatcool=heatcool/SQRT(1.+aconst*aconst)
        flogt=temp/(meanmwt*mhboltz*rho(ib)**gamma1)-entold(ib)-
     &        0.5*dt*(dentdt(ib)+gamma1*heatcool/rho(ib)**gamma1+
     &        dentold(ib)+gamma1*derad(ib)/rho(ib)**gamma1)

        IF(ABS(flogt).GT.10.0*smallnum)
     &     CALL terror(' root error in entcheck ')
C               ------

        RETURN
        END
