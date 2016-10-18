/*
** svn $Id: mct_roms_swan.h 756 2008-09-14 20:18:28Z jcwarner $
***************************************************** Rajesh Kumar   ***
***************************************************** John C. Warner ***
** Copyright (c) 2002-2014 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines couple ROMS/TOMS to the CICE    model using         **
** the Model Coupling Toolkit (MCT).                                  **
**                                                                    **
************************************************************************
*/

      SUBROUTINE initialize_ocn2cic_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize ROMS and CICE models coupling stream.  This is the      !
!  training phase used to constuct MCT parallel interpolators and      !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mct_coupler_params
      USE mod_kinds
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrT, IendT, JstrT, JendT
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Isize, Jsize, MyError
      integer :: i, ic, j, jc, nprocs
      integer :: nRows, nCols, num_sparse_elems
      integer :: cid, cad
      integer, allocatable :: length(:)
      integer, allocatable :: start(:)
      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      character (len=70)    :: nc_name
      character (len=20)    :: to_add
      character (len=120)   :: iostring
      character (len=120)   :: oistring
      real(r8) :: cff
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
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
!
!-----------------------------------------------------------------------
!  Establish MCT communicator.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)
!
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nocn_grids))
        ALLOCATE(AttrVect_G(Nocn_grids))
      END IF
!
!  Initialize MCT coupled model registry.
!
      OCNid=ocnids(ng)
      IF (Nocn_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,myids=ocnids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,OCNid)
      END IF
!
!  Determine start and lengths for domain decomposition.
!
      Jsize=JendR-JstrR+1
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      DO j=JstrR,JendR
        jc=jc+1
        start (jc)=(j)*(Lm(ng)+2)+IstrR+1
        length(jc)=(IendR-IstrR+1)
      END DO
      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapROMS,             &
     &                        start, length, 0, OCN_COMM_WORLD, OCNid)
!
!  Deallocate working arrays.
!
      IF (allocated(start)) THEN
        deallocate (start)
      END IF
      IF (allocated(length)) THEN
        deallocate (length)
      END IF
!
!
!   Initialize attribute vector holding the data code strings from
!   the CICE model.
!
      cad=LEN(iostring)
      DO i=1,cad
        iostring(i:i)=''
      END DO
      cid=1
!
      to_add='AICE'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':freshAI'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':fsaltAI'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':fhocnAI'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':fswthruAI'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':strocnx'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':strocny'
      cad=LEN_TRIM(to_add)
      write(iostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
!  Finalize and remove trailing spaces from the iostring
!  for the rlist.
!
      cad=LEN_TRIM(iostring)
      iostring=iostring(1:cad)
!
!  Initialize attribute vector holding the export data of
!  the CICE model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%cic2ocn_AV,                     &
     &                   rList=TRIM(iostring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%cic2ocn_AV)
!
!  Initialize attribute vector that contain the data strings from
!  the ocean model.
!

      cad=LEN(oistring)
      DO i=1,cad
        oistring(i:i)=''
      END DO
      cid=1
!
      to_add='SST'
      cad=LEN_TRIM(to_add)
      write(oistring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':SSS'
      cad=LEN_TRIM(to_add)
      write(oistring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':FRZMLT'
      cad=LEN_TRIM(to_add)
      write(oistring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
! RAJ
!  u and v needs to cross check VELX and VELY
! RAJ
      to_add=':VELX'
      cad=LEN_TRIM(to_add)
      write(oistring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELY'
      cad=LEN_TRIM(to_add)
      write(oistring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad

!      to_add=':SSH'
      cad=LEN_TRIM(to_add)
      write(oistring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad

!  Finalize and remove trailing spaces from the oistring
!  for the rlist.
!
      cad=LEN_TRIM(oistring)
      oistring=oistring(1:cad)
!
      CALL AttrVect_init(AttrVect_G(ng)%ocn2cic_AV,                     &
     &                   rList=TRIM(oistring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%ocn2cic_AV)
!
      RETURN
      END SUBROUTINE initialize_ocn2cic_coupling

      SUBROUTINE initialize_ocn2cic_routers (tile)
!
!=======================================================================
!                                                                      !
!  Initialize ROMS and CICE models coupling stream.  This is the       !
!  training phase used to constuct MCT parallel interpolators and      !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mct_coupler_params
!
!  Imported variable definitions.
!
      integer, intent(in) :: tile
!
!  Local variable declarations.
!
      integer :: MyError, nprocs
      integer :: ng, ic
!
!-----------------------------------------------------------------------
!  Establish MCT router.
!-----------------------------------------------------------------------
!
      ALLOCATE(Router_C(Nocn_grids,Ncic_grids))
!
!  Initialize routers to the CICE model component.
!
      DO ng=1,Nocn_grids
        DO ic=1,Ncic_grids
          CICid=cicids(ic)
          CALL Router_init (CICid, GlobalSegMap_G(ng)%GSMapROMS,        &
     &                      OCN_COMM_WORLD, Router_C(ng,ic)%ROMStoCICE)
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_ocn2cic_routers

      SUBROUTINE ocn2cic_coupling (ng, ic, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between CICE        !
!  and ROMS  models.                                                   !
!  coded:                                                              !
!                                                                      !
!                                                                      !
!     (...) CICE units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields exported to CICE model:                                      !
!                                                                      !
!     * ----  (), []                          !
!     * ----  (), [ ]                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ic, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2cic_coupling_tile (ng, ic, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif
      RETURN
      END SUBROUTINE ocn2cic_coupling
!
!***********************************************************************
      SUBROUTINE ocn2cic_coupling_tile (ng, ic, tile,                   &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
      USE mct_coupler_params
! CICE
      USE mod_ice
! CICE
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ic, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, MyError, Tag
      integer :: gtype, i, id, ifield, ij, j, k, status

      real(r8) :: add_offset, scale
      real(r8) :: cff, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ubar_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vbar_rho

      real(r8), pointer :: A(:)
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)

      allocate ( A(Asize) )

      A=0.0_r8

!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to seaice (CICE) model.
!-----------------------------------------------------------------------
!
!  SST.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2cic_AV, "SST",    &
     &                           A, Asize)
!
!  SSS.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=OCEAN(ng)%t(i,j,N(ng),nstp(ng),isalt)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2cic_AV, "SSS",    &
     &                           A, Asize)
!
! FRZMLT
!
!jd Assume 5-meter thickness of layer. Here melt potential should be a proper 
!jd weight over a typical mixed (?) / near-ice layer. 

     ij=0
      DO j=JstrR,JendR
       DO i=IstrR,IendR

         if (ICE(ng)%qfraz_accum(i,j).gt.0.0_r8) then
         ICE(ng)%qfraz_accum(i,j)=                                 &
     &        ICE(ng)%qfraz_accum(i,j)/ncouple
             else
       ICE(ng)%qfraz_accum(i,j) = rho0*Cp *                       &
     & min(t_freeze(OCEAN(ng)%t(i,j,N(ng),nstp(ng),isalt),0.0_r8)    &
     &    - OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp), 0.0_r8 )*5.0_r8
         end if

          ij=ij+1
          A(ij)=ICE(ng)%qfraz_accum(i,j)
            end do
           end do

      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2cic_AV,"FRZMLT",  &
     &                           A, Asize)

!jd Reset accumulation array
           ICE(ng)%qfraz_accum(:,:) = 0
!
!
!  U - velocity at RHO-points
! RAJ----needs to check with Kate

     DO j=JstrR,JendR
       DO i=Istr,Iend

         ubar_rho(i,j)=0.5_r8*(OCEAN(ng)%u(i,  j,N(ng),NOUT)+       &
     &      OCEAN(ng)%u(i+1,j,N(ng),NOUT))

        END DO
      END DO

      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          ubar_rho(IstrR,j)=ubar_rho(IstrR+1,j)
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrR,JendR
          ubar_rho(IendR,j)=ubar_rho(IendR-1,j)
        END DO
      END IF

!
!  V-velocity at RHO-points.
!
       DO j=Jstr,Jend
        DO i=IstrR,IendR
        vbar_rho(i,j)=0.5_r8*(OCEAN(ng)%v(i,j,N(ng),NOUT)+        &
      &      OCEAN(ng)%v(i,j+1,N(ng),NOUT))

       END DO
      END DO

      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          vbar_rho(i,JendR)=vbar_rho(i,JendR-1)
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrR,IendR
          vbar_rho(i,JstrR)=vbar_rho(i,JstrR+1)
        END DO
      END IF
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
          A(ij)=ubar_rho(i,j)
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2cic_AV, "VELX",     &
     &                           A, Asize)
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
          A(ij)=vbar_rho(i,j)
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2cic_AV, "VELY",     &
     &                           A, Asize)

!
! RAJ----needs to check with Kate
!

! 
!  SSH
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=OCEAN(ng)%zeta(i,j,KOUT)
        END DO
      END DO

      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2cic_AV, "SSH",    &
     &                           A, Asize)


!
!  Send ocean fields to CICE model.
!
      Tag=ng*100+20+ic

      CALL MCT_isend (AttrVect_G(ng)%ocn2cic_AV,                        &
     &                Router_C(ng,ic)%ROMStoCICE, Tag)

      CALL MCT_waits (Router_C(ng,ic)%ROMStoCICE)

      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'CICE model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,36) ' ** ROMS grid ',ng,                          &
     &                    ' sent data to CICE grid ',ic
 36     FORMAT (a14,i2,a24,i2)
        END IF
      END IF
!
!  Deallocate communication arrays.
!

      deallocate (A)

!
 20   FORMAT (' OCN2CIC_COUPLING - error while sending fields to: ',    &
     &        a, i4)
      RETURN
      END SUBROUTINE ocn2cic_coupling_tile

      SUBROUTINE ocnfcic_coupling (ng, ic, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between cice        !
!  and ocean models.                                                   !
!                                                                      !
!     (...) CICE units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!                                                                      !
!  Fields imported from CICE model:                                    !
!                                                                      !
!     * ---- (), []                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ic, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocnfcic_coupling_tile (ng, ic, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif
      RETURN
      END SUBROUTINE ocnfcic_coupling
!
!***********************************************************************
      SUBROUTINE ocnfcic_coupling_tile (ng, ic, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mct_coupler_params
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
! CICE
      USE mod_ice
! CICE
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ic, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, MyError, Tag
      integer :: gtype, i, id, ifield, ij, j, k, status

      real(r8) :: add_offset, scale
      real(r8) :: cff, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ubar_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vbar_rho

      real(r8), pointer :: A(:)
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)

      allocate ( A(Asize) )

      A=0.0_r8
!
!-----------------------------------------------------------------------
!  Import fields from seaice model (CICE) to ocean model (ROMS).
!-----------------------------------------------------------------------
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      Tag=ng*100+20+ic

!  Receive fields from CICE model.

      CALL MCT_irecv (AttrVect_G(ng)%cic2ocn_AV,                        &
     &                Router_C(ng,ic)%ROMStoCICE, Tag)

!     Wait to make sure the CICE data has arrived.
      CALL MCT_waitr (AttrVect_G(ng)%cic2ocn_AV,                        &
     &                Router_C(ng,ic)%ROMStoCICE)
!
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'CICE model, MyError = ', MyError
        END IF
        exit_flag=2
!       CALL finalize_ocn2cic_coupling
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,38) ' ** ROMS grid ',ng,                          &
     &                    ' recv data from CICE grid ',ic
 38     FORMAT (a14,i2,a26,i2)
        END IF
      END IF
!
! AICE 
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "AICE",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%aice(i,j)=A(ij)
          ELSE
            ICE(ng)%aice(i,j)=ICE(ng)%aice(i,j)+A(ij)
          END IF
        END DO
      END DO
!jd Transfere to u- and v- points for weigthing with stress

           DO j=JstrR,JendR
              DO i=Istr,IendR
                 ICE(ng)%aice_u(i,j)=                                   &
     &                0.5_r8*(ICE(ng)%aice(i-1,j)+ICE(ng)%aice(i,j))
!# ifdef MASKING
!                 ICE(ng)%aice_u(i,j)=ICE(ng)%aice_u(i,j)*umask(i,j)
!# endif
              END DO
           END DO
           DO j=Jstr,JendR
              DO i=IstrR,IendR
                 ICE(ng)%aice_v(i,j)=                                   &
     &                0.5_r8*(ICE(ng)%aice(i,j-1)+ICE(ng)%aice(i,j))
!# ifdef MASKING
!                 ICE(ng)%aice_v(i,j)=ICE(ng)%aice_v(i,j)*vmask(i,j)
!# endif
              END DO
           END DO

!
! freshAI 
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "freshAI",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%freshAI(i,j)=A(ij)
          ELSE
            ICE(ng)%freshAI(i,j)=ICE(ng)%freshAI(i,j)+A(ij)
          END IF
        END DO
      END DO

!
! fsaltAI 
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "fsaltAI",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%fsaltAI(i,j)=A(ij)
          ELSE
            ICE(ng)%fsaltAI(i,j)=ICE(ng)%fsaltAI(i,j)+A(ij)
          END IF
        END DO
      END DO

!
! fhocnAI
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "fhocnAI",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%fhocnAI(i,j)=A(ij)
          ELSE
            ICE(ng)%fhocnAI(i,j)=ICE(ng)%fhocnAI(i,j)+A(ij)
          END IF
        END DO
      END DO

!
! fswthruAI
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "fswthruAI",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%fswthruAI(i,j)=A(ij)
          ELSE
            ICE(ng)%fswthruAI(i,j)=ICE(ng)%fswthruAI(i,j)+A(ij)
          END IF
        END DO
      END DO

!
! strocnx
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "strocnx",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%strx(i,j)=A(ij)
          ELSE
            ICE(ng)%strx(i,j)=ICE(ng)%strx(i,j)+A(ij)
          END IF
        END DO
      END DO

           DO j=JstrR,JendR
              DO i=Istr,IendR
                 ICE(ng)%stru(i,j) = 0.5_r8*                            &
     &                (ICE(ng)%strx(i-1,j) + ICE(ng)%strx(i,j))

!# ifdef MASKING
!                 ICE(ng)%stru(i,j)=ICE(ng)%stru(i,j)*umask(i,j)
!# endif
              END DO
           END DO
!
! strocny
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%cic2ocn_AV, "strocny",   &
     &                           A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          IF (ic.eq.1) THEN
            ICE(ng)%stry(i,j)=A(ij)
          ELSE
            ICE(ng)%stry(i,j)=ICE(ng)%stry(i,j)+A(ij)
          END IF
        END DO
      END DO

          DO j=Jstr,JendR
              DO i=IstrR,IendR
                 ICE(ng)%strv(i,j)=0.5_r8*                              &
     &                (ICE(ng)%stry(i,j-1)+ICE(ng)%stry(i,j))
!# ifdef MASKING
!                 ICE(ng)%strv(i,j)=ICE(ng)%strv(i,j)*vmask(i,j)
!# endif
              END DO
           END DO


!  Apply boundary conditions.


      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!  NEED TO ASK KATE for exchange_u2d_tile and exchange_v2d_tile.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! KATE
!      CALL exchange_r2d_tile (ng, tile,                                 &
!     &                        LBi, UBi, LBj, UBj,                       &
!     &                        ICE(ng)%aice)
! KATE

      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ICE(ng)%aice_u)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ICE(ng)%aice_v)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ICE(ng)%freshAI)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ICE(ng)%fsaltAI)

      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ICE(ng)%fhocnAI)

      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ICE(ng)%fswthruAI)
        CALL exchange_u2d_tile (ng, tile,                              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ICE(ng)%stru)

        CALL exchange_v2d_tile (ng, tile,                               &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ICE(ng)%strv)

      END IF

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
! KATE
!      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    ICE(ng)%aice)
!
! KATE  ---------------

      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ICE(ng)%aice_u, ICE(ng)%aice_v,)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ICE(ng)%freshAI, ICE(ng)%fsaltAI)

      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ICE(ng)%fhocnAI, ICE(ng)%fswthruAI)

      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ICE(ng)%strocnx, ICE(ng)%strocny)

#endif
!
!  Deallocate communication arrays.

      deallocate (A)

!
 10   FORMAT (' OCN2CIC_COUPLING - error while receiving fields from ', &
     &        a, i4)

      RETURN
      END SUBROUTINE ocnfcic_coupling_tile

      SUBROUTINE finalize_ocn2cic_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and seaice models coupling data streams.!
!                                                                       !
!========================================================================
      USE mod_scalars
      USE mct_coupler_params
!
!  Local variable declarations.
!
      integer :: ng, ic, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      deallocate ( cicids )
      deallocate ( ocnids )
      DO ng=1,Nocn_grids
        DO ic=1,Ncic_grids
          CALL Router_clean (Router_C(ng,ic)%ROMStoCICE, MyError)
          CALL AttrVect_clean (AttrVect_G(ng)%ocn2cic_AV, MyError)
        END DO
      END DO

!      CALL GlobalSegMap_clean (GSMapROMS, MyError)

      RETURN

      END SUBROUTINE finalize_ocn2cic_coupling
