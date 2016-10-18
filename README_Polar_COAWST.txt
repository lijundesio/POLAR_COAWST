

For Kate.....

Please check the following..!

1. Master/read_model_inputs.F
 
 	In general this module needs to have access the ice_in (input namelist from CICE) to get the 'dt' and Ncic_grids. 
     Currently we are working on single grid...I think we can hardcode the value of Ncic_grids and read from the CPL (coupling) name list. 
     Also One has to add the same 'dt' value from ice_in to the CPL file each time matching with the ice_in file....

   Please check $/POLAR_COAWST/sample_coupling_inputfile.in


==================================================================================================================================
2. POLAR_COAWST/cice_code/drivers/cice/CICE_RunMod.F90

CICE_Run(coupling_interval)

#ifdef ROMSCOUPLED
      use CICE_MCT, only: CICE_MCT_coupling,TimeInterval

      real(kind=dbl_kind),optional :: coupling_interval

      if (present(coupling_interval)) TimeInterval=coupling_interval
#endif

#ifdef CICE_COUPLING
      IF (MyColor.eq.CICcolor) THEN
        MPI_COMM_ICE=MyCOMM
!        ice_Nmodels=Nmodels
!        ice_CICEid=CICEid
!        ice_OCNid=OCNid

        CALL CICE_INITIALIZE (MyCOMM)
        CALL CICE_RUN(TimeInterval(Iocean,Icice))
!        CALL CICE_RUN
        CALL CICE_FINALIZE

!        CALL CICE_finalize(.TRUE.)
      END IF
#endif

In the above case TimeInterval or the coupling_interval from the namelist (coupling.in) or from other source needs to be passed to CICE model...


==================================================================================================================================

3. for AICE

            ICE(ng)%aice(i,j)=A(ij)

!jd Transfere to u- and v- points for weigthing with stress

                 ICE(ng)%aice_u(i,j)=                                   &
     &                0.5_r8*(ICE(ng)%aice(i-1,j)+ICE(ng)%aice(i,j))


............ ICE(ng)%aice_v(i,j) ...........comp...

==================================================================================================================================

4.   Master/mct_roms_cice.h


 NEED TO ASK KATE for exchange_u2d_tile and exchange_v2d_tile.

==================================================================================================================================

5. home/rajesh/models/COAWST/ROMS/Nonlinear/initial.F

#ifdef CICE_MODEL
      USE CICE_InitMod
#endif
#ifdef CICE_MODEL
! kca - calculate timesteps for restart output
        tspy(ng) = 60*60*24*365/dt(ng)
        tspd(ng) = 60*60*24/dt(ng)
#endif
#ifdef CICE_MODEL
! kca - calculate timesteps for restart output
      dleftinSep = 273-dstart
      Oct = dleftinSep
      Nov = Oct + 31
      Dec = Nov + 30
      Jan = Dec + 31
      Feb = Jan + 31
      Mar = Feb + 28
      Apr = Mar + 31
      May = Apr + 30
      Jun = May + 31
      Jul = Jun + 30
      Aug = Jul + 31
      Sep = Aug + 31
#endif
#ifdef CICE_MODEL
!
!-----------------------------------------------------------------------
! Initialize the ice model
!-----------------------------------------------------------------------
!
      CALL cice_init
#endif

whereas I did......

#ifdef CICE_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled cice model.
!-----------------------------------------------------------------------
!
      allocate(roms_ficoup(Nocn_grids))
      allocate(roms_2icoup(Nocn_grids))
      DO io=1,Nocn_grids
        roms_ficoup(io)=0
        roms_2icoup(io)=0
      END DO
      DO ic=1,Ncic_grids
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocnfcic_coupling (ng, ic, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
      DO ic=1,Ncic_grids
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocn2cic_coupling (ng, ic, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
#endif
!

==================================================================================================================================

6.   /Met_no/roms_code/main3d.F

!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL step3d_t (ng, tile)
# if defined CICE_OCEAN
            CALL frazil_ice_prod (ng,tile)
# endif
          END DO
!$OMP BARRIER
        END DO

==================================================================================================================================

7.  training/POLAR_COAWST/Master/roms_export.F

added 

# ifdef CICE_COUPLING

      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR

      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)

      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
         IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
         IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
         IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
         IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
         JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
         JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
         JendR=BOUNDS(ng)%Jend(tile)+1
      ELSE
         JendR=BOUNDS(ng)%Jend(tile)
      END IF

# else
#  include "set_bounds.h"
# endif

  DO j=JstrR,JendR
          DO i=IstrR,IendR

instead of  DO j=JstrT,JendT
          DO i=IstrT,IendT 

==================================================================================================================================

8 . Met_Roms/apps/common/modified_src/roms-3.6/read_phypar.F

not clear what to do......


==================================================================================================================================

9. /home/rajesh/Desktop/training/POLAR_COAWST/roms_code/Modules/mod_strings.F


#ifdef CICE_OCEAN
     &      'Two-way CICE - ROMS coupling......................',       &
#else
     &      'Two-way Atmosphere-Ocean models coupling .........',       &
#endif


NEED to make it for COAWST....where atm-ocean models coupling is also there....for the time being this is OK since we are only using ROMS-CICE when CICE is defined...!!!

==================================================================================================================================

10. ROMS_KATE/ROMS/Modules/mod_scalars.F

In ROMS-CICE using fake_coupler......below statements are added.....whereas in MetNo....only nudging stuffs are added ...!!!

check with KATE....

#ifdef CICE_MODEL
!    tspy          timesteps per year
!    tspd          timesteps per day
!    dleftinSep    days in September before Oct1 restart
!    Jan,Feb,...,Nov,Dec  days after restart until first of each month
#endif

#ifdef CICE_MODEL
        real(r8), allocatable :: tspy(:)                 ! timesteps
        real(r8), allocatable :: tspd(:)                 ! timesteps
        real(r8) :: dleftinSep = 0.0_r8                  ! days
        real(r8) :: Jan = 0.0_r8                         ! days
        real(r8) :: Feb = 0.0_r8                         ! days
        real(r8) :: Mar = 0.0_r8                         ! days
        real(r8) :: Apr = 0.0_r8                         ! days
        real(r8) :: May = 0.0_r8                         ! days
        real(r8) :: Jun = 0.0_r8                         ! days
        real(r8) :: Jul = 0.0_r8                         ! days
        real(r8) :: Aug = 0.0_r8                         ! days
        real(r8) :: Sep = 0.0_r8                         ! days
        real(r8) :: Oct = 0.0_r8                         ! days
        real(r8) :: Nov = 0.0_r8                         ! days
        real(r8) :: Dec = 0.0_r8                         ! days
        real(r8), allocatable :: rhoice(:)
        real(r8), allocatable :: min_a(:)
#endif

#ifdef CICE_MODEL
      allocate ( tspy(Ngrids) )
      allocate ( tspd(Ngrids) )
      allocate ( rhoice(Ngrids) )
      allocate ( min_a(Ngrids) )
#endif

#ifdef CICE_MODEL
      rhoice = 900.0_r8
      min_a = 1.0e-11_r8
#endif


==================================================================================================================================

11. mod_ice.F 

/training/POLAR_COAWST/roms_code/Modules/mod_ice.F uses the same Met_Roms/apps/common/modified_src/roms-3.6/mod_ice.F....

It does both allocation and initialization in subroutine allocate_ice............which is called in mod_arrays.F

cross check with KATE....

==================================================================================================================================

12. frazil_ice_prod_mod.F 

In Met_Roms/apps/common/modified_src/roms-3.6/...is used insed main3d.F......

Needs to cross check with KATE


==================================================================================================================================

13. ROMS/Nonlinear/bulk_flux.F.....

For LONGWAVE and shortwave...

added according to ROMS-CICE fake_coupler style...

# ifdef CICE_MODEL
     &                     FORCES(ng) % LW_down,                        &
     &                     FORCES(ng) % SW_down,                        &
# endif


NEEDS to cross check this code with KATE....

==================================================================================================================================

14. coupling.dat   ........

==================================================================================================================================



