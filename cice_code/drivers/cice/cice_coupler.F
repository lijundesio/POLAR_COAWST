#if defined ROMSCOUPLED

      MODULE cice_coupler_mod

!
!Modified for CICE-ROMS coupling: cice_coupler.F Author: Rajesh (rajesh.kumar@nyu.edu) 
!
!svn $Id: ocean_coupler.F 755 2008-09-14 19:07:08Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2014 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between        !
!  ROMS/TOMS and other coupled model(s)  via the Model Coupling        !
!  Toolkit (MCT), developed at the Argonne National Laboratory.        !
!                                                                      !
!=======================================================================
!
!
!  Component Model Registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain Decomposition Descriptor DataType and associated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field Storage DataType and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_lsize => lsize
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_copy => copy
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communications scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY: MCT_send => send
      USE m_Transfer, ONLY: MCT_recv => recv
      USE m_Transfer, ONLY: MCT_isend => isend
      USE m_Transfer, ONLY: MCT_irecv => irecv
      USE m_Transfer, ONLY: MCT_waitr => waitrecv
      USE m_Transfer, ONLY: MCT_waits => waitsend
!
! -----------------------------------------------------------------------
!
      USE ice_kinds_mod
      use ice_blocks, only : block, get_block, nx_block, ny_block
      use ice_constants, only: c0, p5, field_loc_center, &
           field_type_scalar, field_loc_NEcorner, &
           field_type_vector
      use ice_domain, only : nblocks, blocks_ice, halo_info
      use ice_domain_size, only : nx_global, ny_global, max_blocks !, block_size_x, block_size_y
      use ice_flux, only: sst, uocn, vocn, zeta, ss_tltx, ss_tlty,&
          sss,frzmlt
      use ice_boundary, only: ice_HaloUpdate
      use ice_fileunits, only: ice_stdout, ice_stderr ! these might be the same

      use ice_accum_shared, only:idaice,idfresh,idfsalt,idfhocn, &
           idfswthru, idstrocnx, idstrocny, accum_time
      use ice_accum_fields, only: accum_i2o_fields, mean_i2o_fields, &
          zero_i2o_fields
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          ice_timer_print, timer_cplrecv, timer_rcvsnd, &
          timer_cplsend,timer_sndrcv,timer_tmp
      USE ice_communicate, ONLY: MPI_COMM_ICE, my_task, master_task

      implicit none
!
      PRIVATE
 
      logical :: initial_call
      real (kind=dbl_kind) :: TimeInterval = 3600.0

      PUBLIC :: initialize_cic_coupling
      PUBLIC :: cic2ocn_coupling
      PUBLIC :: cicfocn_coupling
      PUBLIC :: finalize_cic_coupling
      PUBLIC :: TimeInterval
! Declarations 

! GloabalSegMap variables
      TYPE T_GlobalSegMap_G
          TYPE(GlobalSegMap) :: GSMapCICE      
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:) 
                                                                
! AttrVect variables
      TYPE T_AttrVect_G
        TYPE(AttrVect) :: cic2ocn_AV            ! AttrVect variables
        TYPE(AttrVect) :: ocn2cic_AV
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)

! Router variables
       TYPE T_Router_O
          TYPE(Router)   :: CICEtoROMS           ! Router variables
       END TYPE T_Router_O
       TYPE (T_Router_O), ALLOCATABLE :: Router_O(:,:)

!=======================================================================
!=======================================================================
       INTEGER (int_kind), public :: &
           CICid, OCNid, N_mctmodels
!=======================================================================
!=======================================================================


      CONTAINS

     SUBROUTINE initialize_cic_coupling()
!
!=======================================================================
!                                                                      !
!  Initialize cice and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!

!      implicit none

      include 'mpif.h'
!
!  Local variable declarations.  
!
    integer, pointer :: start(:), length(:)
    integer :: Asize,Istr,Jstr

    integer     :: ilo_glob, j_glob
    integer     :: i, j, iblk, n, gi
    integer     :: lsize,gsize
    integer     :: ier
    integer     :: ic,io
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain

!    character (len=120)   :: iostring
!    character (len=120)   :: oistring

    type(block) :: this_block         ! block information for current block

! RAJ
! Hardcoded
       ic = 1
       io = 1

!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (MPI_COMM_ICE, MyRank, MyError)
      CALL mpi_comm_size (MPI_COMM_ICE, nprocs, MyError)
!
!
! -----------------------------------------------------------------------------------------
!  Initialize MCT coupled model registry.
!
! -----------------------------------------------------------------------------------------
!         ALLOCATE(GlobalSegMap_G(Ncic_grids))  ! if using use mct_couple_params ???
!         ALLOCATE(AttrVect_G(Ncic_grids))
! -----------------------------------------------------------------------------------------
         ALLOCATE(GlobalSegMap_G(ic))
         ALLOCATE(AttrVect_G(ic))

!
!  Initialize MCT coupled model registry.
!
    call ice_timer_start(timer_tmp)

!=======================================================================
!=======================================================================
!    CALL MCTWorld_init (Nmodels, MPI_COMM_WORLD, MPI_COMM_ICE, CICEid)

        CALL MCTWorld_init (N_mctmodels,MPI_COMM_WORLD,                 &
     &                      MPI_COMM_ICE,CICid)

!=======================================================================
!=======================================================================
    WRITE (ice_stdout,*) ' CICE: MCTWorld_init called'

!  Initialize a Global Segment Map for non-haloed transfer of data out
!  of CICE. Determine non-haloed start and length arrays for this
!  processor.
!

! From Met.No.
        n= 0
        
        do iblk = 1, nblocks
        this_block = get_block(blocks_ice(iblk), iblk)
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi
        n = n + (1+jhi-jlo)
        enddo
        lsize = n
        write(ice_stdout, *) ' CICE : lsize  = ', lsize

        allocate (start(lsize))
        allocate (length(lsize))

        do iblk= 1, nblocks
        this_block = get_block(block_ice(iblk), iblk)
       
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi
    
          do j = jlo, jhi
            n = n+1
         ilo_glob = this_block%i_glob(ilo)
         j_glob =  this_block%j_glob(j)
         start(n) = (j_blob-1)*nx_global + ilo_glob
         length(n) = 1+ihi-ilo
          enddo
        enddo
!
!   use grid decomposition to initialize global segmentation map
!
      write (ice_stdout, *) ' CICE : GlobalSegMap_init '

      CALL GlobalSegMap_init (GlobalSegMap_G(ic)%GSMapCICE, start,      &
     &                        length, 0, MPI_COMM_ICE, CICid)
       
      Asize = GlobalSegMap_lsize(GlobalSegMap_G(ic)%GSMapCICE, MPI_COMM_ICE)

      deallocate (start)
      deallocate (length)

!
!  Initialize attribute vector holding the import data code string of
!  the ROMS model.
!

      cad=LEN(ocstring)
      DO i=1,cad
        ocstring(i:i)=''
      END DO
      cid=1
!
      to_add='SST'
      cad=LEN_TRIM(to_add)
      write(ocstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':SSS'
      cad=LEN_TRIM(to_add)
      write(ocstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':FRZMLT'
      cad=LEN_TRIM(to_add)
      write(ocstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELX'
      cad=LEN_TRIM(to_add)
      write(ocstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELY'
      cad=LEN_TRIM(to_add)
      write(ocstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':SSH'
      cad=LEN_TRIM(to_add)
      write(ocstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad

!  Finalize and remove trailing spaces from the ocstring
!  for the rlist.
!
      cad=LEN_TRIM(ocstring)
      ocstring=ocstring(1:cad)
!
!   Initialize import attribute vectors
!
      write (ice_stdout,*) ' CICE : AttrVect_init, Asize = ', Asize
      CALL AttrVect_init(AttrVect_G(ic)%ocn2cic_AV,                     &
     &                   rList=TRIM(ocstring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ic)%ocn2cic_AV)

!
!  Initialize attribute vector holding the export data code string of
!  the CICE model.
!
      cad=LEN(costring)
      DO i=1,cad
        costring(i:i)=''
      END DO
      cid=1
!
! 
      to_add='AICE'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':freshAI'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':fsaltAI'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':fhocnAI'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':fswthruAI'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':strocnx'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':strocny'
      cad=LEN_TRIM(to_add)
      write(costring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
!  Finalize and remove trailing spaces from the costring
!  for the rlist.
!
      cad=LEN_TRIM(costring)
      costring=costring(1:cad)
!
!
!   Initialize export attribute vectors
!
      CALL AttrVect_init (AttrVect_G(ic)%cic2ocn_AV,                    &
     &                    rList=TRIM(costring),lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ic)%cic2ocn_AV)
!

!
!  Initialize a router to the ocean model component.
!
      write (ice_stdout,*) ' CICE : Router_init '

!  RAJ ...not dealing with nested grid for the time being....!  
!      allocate (Router_O(Ncic_grids,Nocn_grids))
      allocate (Router_O(ic,io))
!  RAJ  

      call Router_init (OCNid, GlobalSegMap_G(ic)%GSMapCICE,       &
                        MPI_COMM_ICE, Router_O(ic,io)%CICEtoROMS)

      write (ice_stdout,*) ' CICE : Router_init.  DONE '
! -------------------------
      initial_call = .true.

      call ice_timer_stop (timer_tmp)
      call ice_timer_print (timer_tmp, .false.)


      RETURN
      END SUBROUTINE initialize_cic_coupling


      SUBROUTINE cic2ocn_coupling ()
!
!=======================================================================
!                                                                      
!  This subroutine reads and writes the coupling data streams between  
!  cice and ocean models. Currently, the following data streams are    
!  processed:                                                          
!                                                                      
!  Fields exported to the OCEAN model:field (unit)                                 
!                                                                      
!     * AICE ()                                   
!     * Fresh AI ()                                        
!     * Fsalt AI ()                          
!     * Fhocn AI ()                                         
!     * Fswthru AI ()                      
!     * Strocnx and Strocny                       
!                                                                      
!  Fields imported from the OCEAN Model:                               
!                                                                      
!     * SST ()
!     * SSS ()
!     * FRZMLT ()
!     * VELX and VELY (m/s)      
!     * SSH (m) 
!=======================================================================
!
      USE ice_grid, only : HTN, HTE, dxu, dyu, dxt, dyt, &
         t2ugrid_vector
      use ice_calendar, only: dt, time, write_ic ,istep, istep1
 
      implicit none

    real(kind=dbl_kind), pointer :: avdata(:)
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer     :: i,j,Csize,iblk,n
    integer     :: ic,io,Tag

    integer     :: MyRank,MyError,NPROCS

! RAJ
! Hardcoded for single grid 
       ic = 1
       io = 1
!
!-----------------------------------------------------------------------
!  Get communicator local rank and size
!-----------------------------------------------------------------------
!
! RAJ+++++++++++++++++++++++++++
      CALL MPI_COMM_RANK (MPI_COMM_ICE, MyRank, MyError)
      CALL MPI_COMM_SIZE (MPI_COMM_ICE, NPROCS, MyError)
! RAJ+++++++++++++++++++++++++++
!
!        ***********************************
!          Send CICE fields to ROMS 
!        ***********************************
!

!   update of accumulated time since last coupling (accum_time) and accumulation
!   of flux fields are done from RunMod (calling functions in ice_accum_fields).

       IF (accum_time >= TimeInterval .or. initial_call) THEN
 
         call ice_timer_start(timer_sndrcv)
       IF (my_task == master_task) THEN
          write(ice_stdout,*) '*********************************************'

          write(ice_stdout,*) 'CICE - Ocean: coupling routine called from CICE'
          write(ice_stdout,*) 'time = ', time
          write(ice_stdout,*) 'dt = ', dt
          write(ice_stdout,*) 'istep ', istep
          write(ice_stdout,*) 'istep1 ', istep1
          write(ice_stdout,*) '*********************************************'
          call flush(ice_stdout)
       END IF

!  Get the number of grid point on this processor.
!
      Csize=GlobalSegMap_lsize(GlobalSegMap_G(ic)%GSMapCICE,          &
     &                           MPI_COMM_ICE)
!
       allocate(avdata(Csize))
       avdata=0.0

       call mean_i2o_fields()
!
!-----------------------------------------------------------------------
!  Export fields from cice (CICE) to ocean (ROMS) model.
!
!     AICE        sea ice concentration.
!
!
      field=accum_i2o_fields(:,:,idaice,:)    ! AICE

      call field2avec(field,avdata)       

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV, "AICE",     &
     &                           avdata)
!
!      Fresh_ai

      fied=accum_i2o_fields(:,:,idfresh,:)
      call field2avec(field,avdata)

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV, "freshAI",  &
     &                           avdata)
!
!      FsaltAI

      fied=accum_i2o_fields(:,:,idfsalt,:)
      call field2avec(field,avdata)

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV, "fsaltAI",  &
     &                           avdata)
!
!      FhocnAI

      fied=accum_i2o_fields(:,:,idfhocn,:)
      call field2avec(field,avdata)

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV, "fhocnAI",  &
     &                           avdata)
!
!      FswthruAI
      fied=accum_i2o_fields(:,:,idfswthru,:)
      call field2avec(field,avdata)

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV,"fswthruAI", &
     &                           avdata)
!

! Export stress vector. Allready converted to T-cell and scaled with aice
! Change of sign here as the stress on the ocean acts in opposite
! directon as the stress on the ice.

!      strocnx
      fied=-accum_i2o_fields(:,:,idstrocnx,:)
      call field2avec(field,avdata)

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV,"strocnx",   &
     &                           avdata)
!

!      strocnx
      fied=-accum_i2o_fields(:,:,idstrocny,:)
      call field2avec(field,avdata)

      CALL AttrVect_importRAttr (AttrVect_G(ic)%cic2ocn_AV,"strocny",   &
     &                           avdata)
!

      call ice_timer_stop(timer_sndrcv)

      Tag=io*100+20+ic

! Transfere data to ocean
       call ice_timer_start(timer_cplsend)
!       CALL MCT_Send(cice2ocn_AV, CICEtoROMS)
!RAJ
         CALL MCT_isend (AttrVect_G(ic)%cic2ocn_AV,                    &
     &                   Router_O(ic,io)%CICEtoROMS, Tag)

          CALL MCT_waits (Router_O(ic,io)%CICEtoROMS)
!RAJ

       call ice_timer_stop(timer_cplsend)
       if (initial_call) then
          call ice_timer_print(timer_cplsend,.false.)
       endif


          IF (my_task == master_task) THEN
            WRITE (ice_stdout,*)'CICE : sent data to ROMS grid '
          END IF

          IF (MyError.ne.0) THEN
       WRITE (*,*)'coupling send fail cice_coupler, Error= ', MyError
            CALL FINALIZE_CIC_COUPLING(ic)
          END IF

      deallocate (avdata)
!
!     end if for---- IF (accum_time >= TimeInterval .or. initial_call) THEN
        END IF 
       initial_call=.false.

!        ***********************************

contains
    subroutine field2avec(field, avec)

    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer     :: i,j,Asize,iblk,n

      real(kind=dbl_kind), pointer :: avec(:)
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
           intent(in):: field
      n = 0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
            do i = ilo, ihi
               n = n+1
               avec(n) = field(i,j,iblk)
            enddo
         enddo
      enddo
    end subroutine field2avec

       RETURN
      END SUBROUTINE cic2ocn_coupling


         SUBROUTINE cicfocn_coupling ()
!
!=======================================================================
!                                                                      
!  This subroutine reads and writes the coupling data streams between  
!  cice and ocean models. Currently, the following data streams are    
!  processed:                                                          
!                                                                      
!  Fields exported to the OCEAN model:field (unit)                                 
!                                                                      
!     * AICE ()                                   
!     * Fresh AI ()                                        
!     * Fsalt AI ()                          
!     * Fhocn AI ()                                         
!     * Fswthru AI ()                      
!     * Strocnx and Strocny                       
!                                                                      
!  Fields imported from the OCEAN Model:                               
!                                                                      
!     * SST ()
!     * SSS ()
!     * FRZMLT ()
!     * VELX and VELY (m/s)      
!     * SSH (m) 
!=======================================================================

      USE ice_grid, only : HTN, HTE, dxu, dyu, dxt, dyt, &
         t2ugrid_vector
      use ice_calendar, only: dt, time, write_ic ,istep, istep1

      implicit none

    real(kind=dbl_kind), pointer :: avdata(:)
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer     :: i,j,Csize,iblk,n
    integer     :: ic,io,Tag

! RAJ
! Hardcoded for single grid 
       ic = 1
       io = 1

!-----------------------------------------------------------------------
!  Get communicator local rank and size
!-----------------------------------------------------------------------

      CALL MPI_COMM_RANK (MPI_COMM_ICE, MyRank, MyError)
      CALL MPI_COMM_SIZE (MPI_COMM_ICE, NPROCS, MyError)
!
       IF (accum_time >= TimeInterval .or. initial_call) THEN
!  Get the number of grid point on this processor.
!
      Csize=GlobalSegMap_lsize(GlobalSegMap_G(ic)%GSMapCICE,          &
     &                           MPI_COMM_ICE)

!  Allocate attribute vector array used to export/import data.
!
       allocate (avdata(Csize))
       avdata=0.0

! RAJ -- needs to check reqd or not....used in cic2ocn...!
        call mean_i20_fields()
! RAJ -- needs to check reqd or not....used in cic2ocn...!

!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) model CICE  model
!-----------------------------------------------------------------------
!  Receive data from ocean
!
          Tag=io*100+20+ic
 
       call ice_timer_start(timer_cplrecv)

          CALL MCT_irecv (AttrVect_G(ic)%ocn2cic_AV,                    &
     &                   Router_O(ic,io)%CICEtoROMS, Tag)

          CALL MCT_waitr (AttrVect_G(ic)%ocn2cic_AV,                    &
     &                 Router_O(ic,io)%CICEtoROMS)

       call ice_timer_stop(timer_cplrecv)


       call ice_timer_start(timer_rcvsnd)

          IF (my_task==master_task) THEN
       WRITE (ice_stdout,*)'CICE grid  recv data from ROMS grid '
          END IF

          IF (MyError.ne.0) THEN
       WRITE (*,*)'cice_coupler err-No data from ROMS,Error= ',MyError
            CALL FINALIZE_CIC_COUPLING(ic)
          END IF
!
!     SST

      CALL AttrVect_exportRAttr (AttrVect_G(ic)%ocn2cic_AV,             &
     &                           'SST', avdata)

       write(ice_stdout,*) 'CICE rank ',my_task,  &
            ' setting the sst field (max/min): ', &
            maxval(avdata), ' ', minval(avdata)
       call avec2field(avdata,sst)
       call ice_HaloUpdate (sst, halo_info, &
            field_loc_center, field_type_scalar)

! Salinity
       CALL AttrVect_exportRAttr(AttrVect_G(ic)%ocn2cice_AV,      &
     &            'SSS', avdata)

       write(ice_stdout,*) 'CICE rank ',my_task,  &
            ' setting the sss field (max/min): ', &
            maxval(avdata), ' ', minval(avdata)

       call avec2field(avdata,sss)

       if (minval(avdata) < c0) then
          write(ice_stdout,*) 'CICE rank ',my_task,  &
               ' correcting invalid sss ', minval(avdata)

          do iblk = 1, nblocks
             this_block = get_block(blocks_ice(iblk),iblk)
             ilo = this_block%ilo
             ihi = this_block%ihi
             jlo = this_block%jlo
             jhi = this_block%jhi
             do j = jlo, jhi
                do i = ilo, ihi
                   sss(i,j,iblk)=max(sss(i,j,iblk),c0)
                end do
             end do
          end do
       endif
       call ice_HaloUpdate (sss, halo_info, &
            field_loc_center, field_type_scalar)

! Melt freeze potential
       CALL AttrVect_exportRAttr(AttrVect_G(ic)%ocn2cice_AV,     &
     &          'FRZMLT', avdata)

       write(ice_stdout,*) 'CICE rank ',my_task,  &
            ' setting the frzmlt field (max/min): ', &
            maxval(avdata), ' ', minval(avdata)

       call avec2field(avdata,frzmlt)
       call ice_HaloUpdate (frzmlt, halo_info, &
            field_loc_center, field_type_scalar)

!
! recieve ocean currents and interpolate to B grid
!
! Uocn
      
       CALL AttrVect_exportRAttr(AttrVect_G(ic)%ocn2cice_AV,   &
     &            'VELX', avdata)

       write(ice_stdout,*) 'CICE rank ', my_task,    &
            ' setting the VELX (Uocn) field(max/min): ',&
            maxval(avdata), ' ', minval(avdata)

       call avec2field(avdata,uocn)

       call t2ugrid_vector(uocn)

       call ice_HaloUpdate (uocn, halo_info, &
            field_loc_NEcorner, field_type_vector)

! Vocn

       CALL AttrVect_exportRAttr(AttrVect_G(ic)%ocn2cice_AV,   &
     &            'VELY', avdata)

       write(ice_stdout,*) 'CICE rank ', my_task, &
            ' setting the VELY (vocn) field(max/min): ', &
            maxval(avdata), ' ', minval(avdata)

       call avec2field(avdata,vocn)
       call ice_HaloUpdate (vocn, halo_info, &
            field_loc_center, field_type_scalar)

       call t2ugrid_vector(vocn)

       call ice_HaloUpdate (vocn, halo_info, &
            field_loc_NEcorner, field_type_vector)
!
! SSH
!
       CALL AttrVect_exportRAttr(AttrVect_G(ic)%ocn2cice_AV,    &
     &            'SSH', avdata)

       write(ice_stdout,*) 'CICE rank ', my_task, &
            ' setting the SSH field(max/min): ', &
            maxval(avdata), ' ', minval(avdata)

       call avec2field(avdata,zeta)
       call ice_HaloUpdate (zeta, halo_info, &
            field_loc_center, field_type_scalar)

       do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          do j = jlo, jhi
             do i = ilo, ihi
                ss_tltx(i,j,iblk) =                                  &
                     p5*((zeta(i+1,j,iblk) + zeta(i+1,j+1,iblk)) -   &
                     (zeta(i,j,iblk) + zeta(i,j+1,iblk))) /dxu(i,j,iblk)
                ss_tlty(i,j,iblk) =                                  &
                     p5*((zeta(i,j+1,iblk) + zeta(i+1,j+1,iblk)) -   &
                     (zeta(i,j,iblk) + zeta(i+1,j,iblk))) /dyu(i,j,iblk)
             enddo
          enddo
       enddo

       call ice_HaloUpdate (ss_tltx, halo_info, &
            field_loc_NEcorner, field_type_vector)
       call ice_HaloUpdate (ss_tlty, halo_info, &
            field_loc_NEcorner, field_type_vector)


       call zero_i2o_fields ! also accum_time is zeroed

       deallocate(avdata)
       call ice_timer_stop(timer_rcvsnd)

    END IF
    initial_call=.false.

!        ***********************************

contains
    subroutine avec2field(avec, field)

    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    integer     :: i,j,iblk,n

      real(kind=dbl_kind), pointer :: avec(:)
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
           intent(out) :: field
      n = 0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
            do i = ilo, ihi
               n = n+1
               field(i,j,iblk)=avec(n)
            enddo
         enddo
      enddo
    end subroutine avec2field

      RETURN
      END SUBROUTINE cicfocn_coupling

!--------------------------------------------------------------
    subroutine finalize_cic_coupling

    implicit none
    integer :: err, ic, io

!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (Router_O(ic,io)%CICEtoROMS, err)
      CALL AttrVect_clean (Router_O(ic,io)%ocn2cic_AV, err)
      CALL AttrVect_clean (Router_O(ic,io)%cic2ocn_AV, err)
      CALL GlobalSegMap_clean (GlobalSegMap_G(ic)%GSMapCICE, err)

  end subroutine finalize_cic_coupling

END MODULE cice_coupler_mod

#endif
