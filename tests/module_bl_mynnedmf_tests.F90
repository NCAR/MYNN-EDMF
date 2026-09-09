! Module for MYNN EDMF scheme tests
module module_bl_mynnedmf_tests
    use module_bl_mynnedmf_driver
    use netcdf
    use, intrinsic :: ieee_exceptions
    use, intrinsic :: ieee_arithmetic
    ! public
    !=================================================================================================================    
    implicit none
    logical, dimension(5) :: halt_flags 
    logical :: flag_iter,bl_mynn_tkeadvect,cycling
    integer :: initflag
    real,dimension(1) :: pattern_spp_pbl
    integer :: bl_mynn_output,spp_pbl
    logical :: &
      flag_qc,               &     ! if true,the physics package includes the cloud liquid water mixing ratio.
      flag_qi,               &     ! if true,the physics package includes the cloud ice mixing ratio.
      flag_qs,               &     ! if true,the physics package includes the snow mixing ratio.
      flag_qnc,              &     ! if true,the physics package includes the cloud liquid water number concentration.
      flag_qni,              &     ! if true,the physics package includes the cloud ice number concentration.
      flag_qnifa,            &     ! if true,the physics package includes the "ice-friendly" aerosol number concentration.
      flag_qnwfa,            &     ! if true,the physics package includes the "water-friendly" aerosol number concentration.
      flag_qnbca                   ! if true,the physics package includes the number concentration of black carbon.
    logical, parameter :: flag_oz = .false.
    
    contains

    subroutine init_mynn_edmf_flags()
       
       write(*,*) '--- calling  init_mynn_edmf_flags() ---'      
       
       cycling=.false.
       initflag=1
       spp_pbl=1
       pattern_spp_pbl=0.0
       flag_iter = .true.
       flag_qc=.true.
       flag_qi=.true.
       flag_qs=.true.
       flag_qnc=.false.
       flag_qni=.false.
       flag_qnifa=.false.
       flag_qnwfa=.false.
       flag_qnbca=.false.
       bl_mynn_tkeadvect=.true.
       bl_mynn_output=1                              

    end subroutine init_mynn_edmf_flags

    !=================================================================================================================    
    subroutine mynnedmf_test(case,bl_mynn_closure,bl_mynn_cloudpdf,bl_mynn_mixlength,      &
        bl_mynn_edmf,bl_mynn_edmf_dd,bl_mynn_edmf_mom,bl_mynn_edmf_tke,bl_mynn_cloudmix,   &
        bl_mynn_mixqt, bl_mynn_mixscalars, bl_mynn_mixaerosols,bl_mynn_mixnumcon,          &
        bl_mynn_ess, tke_budget, bl_mynn_diags, mix_chem, enh_mix, restart_in,             &
        t_start_in, t_end_in, u, v, th, qv, qc, qi,                                        &
        rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten,                          &
        qc_bl, qi_bl, cldfra_bl, el_pbl, qke, qsq, tsq, cov,                               &
        sh, sm, qke_adv, pblh)       
        
        implicit none
        
        character(len=*),intent(in) :: case
        integer :: ncid, varid
        integer :: dimid_time, dimid_z
        integer :: nt, nz
        integer :: t, t_start, t_end, k
        integer :: status
        integer :: ims,ime,kms,kme,jms,jme
        integer :: ids,ide,kds,kde,jds,jde
        integer :: its,ite,kts,kte,jts,jte
        ! integer :: ndims,dimids(10)
 
        logical :: bl_mynn_tkeadvect, cycling, restart1
        integer :: bl_mynn_cloudpdf,                            &
                 bl_mynn_mixlength,                             &
                 bl_mynn_edmf,                                  &
                 bl_mynn_edmf_dd,                               &
                 bl_mynn_edmf_mom,                              &
                 bl_mynn_edmf_tke,                              &
                 bl_mynn_cloudmix,                              &
                 bl_mynn_mixqt,                                 &
                 bl_mynn_output,                                &
                 bl_mynn_mixscalars,                            &
                 bl_mynn_mixaerosols,                           &
                 bl_mynn_mixnumcon,                             &
                 bl_mynn_ess,                                   &
                 spp_pbl,                                       &
                 tke_budget,                                    &
                 bl_mynn_diags
        real ::  bl_mynn_closure
        real :: delt,dxc

        LOGICAL :: ALLOWED_TO_READ
        INTEGER :: P_QC,P_QI,PARAM_FIRST_SCALAR

        character(len=19), allocatable :: time(:)

        ! --- optional interface for restart ---
        integer, intent(in), optional :: t_start_in, t_end_in
        logical, intent(in), optional :: restart_in

        real, dimension(:,:,:), allocatable, intent(inout),  optional :: u,v,th,qv,qc,qi
        real, dimension(:,:,:), allocatable, intent(inout),  optional :: rublten,rvblten,rthblten,  &
                                                                        rqvblten,rqcblten,rqiblten
        real, dimension(:,:,:), allocatable, intent(inout),  optional :: qc_bl,qi_bl,cldfra_bl
        real, dimension(:,:,:), allocatable, intent(inout),  optional :: el_pbl,qke,qsq,tsq,cov,    &
                                                                        sh,sm,qke_adv
        real, dimension(:,:), allocatable, intent(inout),  optional :: pblh
        logical :: is_warm_start

        ! 2D arrays
        real, allocatable :: xland(:,:), ps(:,:), ts(:,:), qsfc(:,:), ust(:,:), ch(:,:),            &
                hfx(:,:), qfx(:,:), wspd(:,:), znt(:,:), uoce(:,:), voce(:,:), dx2d(:,:)
        
        ! output 2D arrays
        real, allocatable :: excess_h(:,:), excess_q(:,:), maxmf(:,:),maxwidth(:,:),                &
             ztop_plume(:,:), maxwidth_dd(:,:), maxmf_dd(:,:), ent_eff(:,:),                        &
             maxtkeprod(:,:), cldtop_cooling(:,:)
        integer, allocatable :: kpbl(:,:)
        real, allocatable :: pblh_loc(:,:)
        real, allocatable :: lwp(:,:), iwp(:,:), swp(:,:), wspd10(:,:), wspd80(:,:), wspd160(:,:),  &
             cldceil(:,:), maxcldfra(:,:), maxcldfra_pbl(:,:)
        
        ! 3D arrays
        real, allocatable :: u_loc(:,:,:),v_loc(:,:,:), w_loc(:,:,:), th_loc(:,:,:), t3d_loc(:,:,:),& 
                p_loc(:,:,:), exner_loc(:,:,:), rho_loc(:,:,:), qv_loc(:,:,:), qc_loc(:,:,:),       &
                qi_loc(:,:,:), dz_loc(:,:,:), exch_h_loc(:,:,:), exch_m_loc(:,:,:),                 &
                pattern_spp_pbl(:,:,:), pint(:,:,:)
        real, allocatable ::  rthraten_loc(:,:,:), rublten_loc(:,:,:), rvblten_loc(:,:,:), rthblten_loc(:,:,:)        
        
        !optional CHEM arrays
        real, allocatable ::  chem3d(:,:,:,:), settle3d(:,:,:,:), vd3d(:,:,:)
        real, allocatable ::  frp_mean(:,:)  , emis_ant_no(:,:)
        
        !optional and output 3D arrays
        real, allocatable :: qc_bl_loc(:,:,:), qi_bl_loc(:,:,:), cldfra_bl_loc(:,:,:)
        real, allocatable :: qke_loc(:,:,:), qke_adv_loc(:,:,:), el_pbl_loc(:,:,:), sh3d_loc(:,:,:),            &
                sm3d_loc(:,:,:), tsq_loc(:,:,:), qsq_loc(:,:,:), cov_loc(:,:,:)
        real, allocatable :: qnbca_loc(:,:,:),qnc_loc(:,:,:),qni_loc(:,:,:),qnifa_loc(:,:,:),qnwfa_loc(:,:,:),  & 
                qs_loc(:,:,:),qshear_loc(:,:,:),qwt_loc(:,:,:),qBUOY_loc(:,:,:),qDISS_loc(:,:,:)
        real, allocatable :: rqcblten_loc(:,:,:),rqiblten_loc(:,:,:),rqnbcablten_loc(:,:,:),                    &
                rqniblten_loc(:,:,:), rqnifablten_loc(:,:,:), rqnwfablten_loc(:,:,:), rqsblten_loc(:,:,:),      &
                rqvblten_loc(:,:,:), sub_sqv3d_loc(:,:,:), rqncblten_loc(:,:,:), sub_thl3d_loc(:,:,:)

        ! output 3D arrays
        real, allocatable :: det_sqv3d_loc(:,:,:),dqke_loc(:,:,:),edmf_a_loc(:,:,:),edmf_ent_loc(:,:,:),        &
                edmf_qc_loc(:,:,:), edmf_qt_loc(:,:,:),edmf_thl_loc(:,:,:),edmf_w_loc(:,:,:),det_thl3d_loc(:,:,:)

        !smoke/dust configurations
        logical, intent(in), optional:: mix_chem
        logical, intent(in), optional:: enh_mix
        integer,parameter:: nchem=1,ndvel=1

        !ccpp obligation
        character::errmsg
        integer::errflg

        ! Open NetCDF file
        print*,'Case: ',trim(case)
        ! Save current halting mode
        call ieee_get_halting_mode(ieee_all, halt_flags)
        
        ! Disable FPE traps for NetCDF operations
        call ieee_set_halting_mode(ieee_all, .false.)
  
        status = nf90_open('./data/input_'//trim(case)//'.nc', NF90_NOWRITE, ncid)
        print*,'status',status
        if (status /= nf90_noerr) then
            print *, "Error opening file: ./data/input_", trim(case),'.nc'
            print *, trim(nf90_strerror(status))
            stop
        endif

        ! Restore original halting mode
        call ieee_set_halting_mode(ieee_all, halt_flags)      
        
        ! Get dimensions
        status = nf90_inq_dimid(ncid, "Time", dimid_time)
        status = nf90_inquire_dimension(ncid, dimid_time, len=nt)
        
        status = nf90_inq_dimid(ncid, "bottom_top", dimid_z)
        status = nf90_inquire_dimension(ncid, dimid_z, len=nz)
        
        print *, "Dimensions: nz=", nz, " nt=", nt

        t_start = 2
        if (present(t_start_in)) t_start = t_start_in
        t_end = nt
        if (present(t_end_in)) t_end = t_end_in

        restart1 = .false.
        if (present(restart_in)) restart1 = restart_in
        is_warm_start = restart1 .and. present(u)

        print *, 'is warm start ', is_warm_start
        print *, 'restart ', restart1

        ims = 1
        ime = 1
        jms = 1
        jme = 1
        kms = 1
        kme = nz+1

        its = 1
        ite = 1
        jts = 1
        jte = 1
        kts = 1
        kte = nz

        ids = 1
        ide = 2
        jds = 1
        jde = 2
        kds = 1
        kde = nz
        delt = 18.
        dxc = 3000.

        ! allocate(time(nt))
        allocate(character(len=19) :: time(nt))        
        
        ! allocate 2D arrays
        allocate(xland(ims:ime, jms:jme))
        allocate(ps(ims:ime, jms:jme))
        allocate(ts(ims:ime, jms:jme))
        allocate(qsfc(ims:ime, jms:jme))
        allocate(ust(ims:ime, jms:jme))
        allocate(ch(ims:ime, jms:jme))
        allocate(hfx(ims:ime, jms:jme))
        allocate(qfx(ims:ime, jms:jme))
        allocate(wspd(ims:ime, jms:jme))
        allocate(znt(ims:ime, jms:jme))
        allocate(uoce(ims:ime, jms:jme))
        allocate(voce(ims:ime, jms:jme))
        allocate(dx2d(ims:ime, jms:jme))

        allocate(kpbl(ims:ime, jms:jme))
        allocate(maxmf(ims:ime, jms:jme))
        allocate(maxwidth(ims:ime, jms:jme))
        allocate(pblh_loc(ims:ime, jms:jme))
        allocate(ztop_plume(ims:ime, jms:jme))
        allocate(excess_h(ims:ime, jms:jme))
        allocate(excess_q(ims:ime, jms:jme))
        allocate(maxwidth_dd(ims:ime, jms:jme))
        allocate(maxmf_dd(ims:ime, jms:jme))
        allocate(ent_eff(ims:ime, jms:jme))
        allocate(maxtkeprod(ims:ime, jms:jme))
        allocate(cldtop_cooling(ims:ime, jms:jme))

        ! allocate 2D diagnostic fields
        allocate(lwp(ims:ime, jms:jme))
        allocate(iwp(ims:ime, jms:jme))
        allocate(swp(ims:ime, jms:jme))
        allocate(wspd10(ims:ime, jms:jme))
        allocate(wspd80(ims:ime, jms:jme))
        allocate(wspd160(ims:ime, jms:jme))
        allocate(cldceil(ims:ime, jms:jme))
        allocate(maxcldfra(ims:ime, jms:jme))
        allocate(maxcldfra_pbl(ims:ime, jms:jme))

        ! allocate 3D arrays
        allocate(qBUOY_loc(ims:ime, kms:kme, jms:jme))
        allocate(qDISS_loc(ims:ime, kms:kme, jms:jme))
        allocate(u_loc(ims:ime, kms:kme, jms:jme))
        allocate(v_loc(ims:ime, kms:kme, jms:jme))
        allocate(w_loc(ims:ime, kms:kme, jms:jme))
        allocate(th_loc(ims:ime, kms:kme, jms:jme))
        allocate(t3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(p_loc(ims:ime, kms:kme, jms:jme))
        allocate(pint(ims:ime, kms:kme, jms:jme))
        allocate(exner_loc(ims:ime, kms:kme, jms:jme))
        allocate(rho_loc(ims:ime, kms:kme, jms:jme))
        allocate(qv_loc(ims:ime, kms:kme, jms:jme))
        allocate(qc_loc(ims:ime, kms:kme, jms:jme))
        allocate(qi_loc(ims:ime, kms:kme, jms:jme))
        allocate(dz_loc(ims:ime, kms:kme, jms:jme))
        allocate(rthraten_loc(ims:ime, kms:kme, jms:jme))
        allocate(rublten_loc(ims:ime, kms:kme, jms:jme))
        allocate(rvblten_loc(ims:ime, kms:kme, jms:jme))
        allocate(rthblten_loc(ims:ime, kms:kme, jms:jme))
        allocate(exch_h_loc(ims:ime, kms:kme, jms:jme))
        allocate(exch_m_loc(ims:ime, kms:kme, jms:jme))

        allocate(cov_loc(ims:ime, kms:kme, jms:jme))
        allocate(det_thl3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(det_sqv3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(dqke_loc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_a_loc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_ent_loc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_qc_loc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_qt_loc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_thl_loc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_w_loc(ims:ime, kms:kme, jms:jme))        
        
        allocate(qnc_loc(ims:ime, kms:kme, jms:jme))  
        allocate(qni_loc(ims:ime, kms:kme, jms:jme))  
        allocate(qnwfa_loc(ims:ime, kms:kme, jms:jme))  
        allocate(qnifa_loc(ims:ime, kms:kme, jms:jme))  
        allocate(qs_loc(ims:ime, kms:kme, jms:jme))
        allocate(qshear_loc(ims:ime, kms:kme, jms:jme)) 
        allocate(qwt_loc(ims:ime, kms:kme, jms:jme))  
        allocate(qke_loc(ims:ime, kms:kme, jms:jme))
        allocate(qke_adv_loc(ims:ime, kms:kme, jms:jme))
        allocate(tsq_loc(ims:ime, kms:kme, jms:jme))
        allocate(qsq_loc(ims:ime, kms:kme, jms:jme))
        allocate(el_pbl_loc(ims:ime, kms:kme, jms:jme))
        allocate(sh3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(sm3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(qc_bl_loc(ims:ime, kms:kme, jms:jme))
        allocate(qi_bl_loc(ims:ime, kms:kme, jms:jme))
        allocate(cldfra_bl_loc(ims:ime, kms:kme, jms:jme))
        allocate(sub_thl3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(sub_sqv3d_loc(ims:ime, kms:kme, jms:jme))
        allocate(RQVBLTEN_loc(ims:ime, kms:kme, jms:jme))
        allocate(RQCBLTEN_loc(ims:ime, kms:kme, jms:jme))
        allocate(RQIBLTEN_loc(ims:ime, kms:kme, jms:jme)) 

        ! allocate psudo smoke/chem test arrays
        allocate(chem3d(ims:ime, kms:kme, jms:jme, 1:nchem))
        allocate(settle3d(ims:ime, kms:kme, jms:jme, 1:nchem))
        allocate(vd3d(ims:ime, kms:kme, jms:jme))
        allocate(frp_mean(ims:ime, jms:jme))
        allocate(emis_ant_no(ims:ime, jms:jme))


        ! Read time variable
        status = nf90_inq_varid(ncid, "Times", varid)
        status = nf90_get_var(ncid, varid, time)
        
        call init_mynn_edmf_flags()
     
        call mynnedmf_init (                               &
           &  RUBLTEN_loc,RVBLTEN_loc,RTHBLTEN_loc,        &
           &  RQVBLTEN_loc,RQCBLTEN_loc,                   &
           &  RQIBLTEN_loc,QKE_loc,                        &
           &  restart1,.true.,                             &
           &  P_QC,P_QI,PARAM_FIRST_SCALAR,                &
           &  IDS,IDE,JDS,JDE,KDS,KDE,                     &
           &  IMS,IME,JMS,JME,KMS,KME,                     &
           &  ITS,ITE,JTS,JTE,KTS,KTE                      )

        ! Copy input values if this is a warm restart
        if (is_warm_start) then
            if (present(u)) then; if (allocated(u)) u_loc(1,:,1) = u(1,:,1); end if
            if (present(v)) then; if (allocated(v)) v_loc(1,:,1) = v(1,:,1); end if
            if (present(th)) then; if (allocated(th)) th_loc(1,:,1) = th(1,:,1); end if
            if (present(qv)) then; if (allocated(qv)) qv_loc(1,:,1) = qv(1,:,1); end if
            if (present(qc)) then; if (allocated(qc)) qc_loc(1,:,1) = qc(1,:,1); end if
            if (present(qi)) then; if (allocated(qi)) qi_loc(1,:,1) = qi(1,:,1); end if

            if (present(rublten)) then; if (allocated(rublten)) rublten_loc(1,:,1) = rublten(1,:,1); end if
            if (present(rvblten)) then; if (allocated(rvblten)) rvblten_loc(1,:,1) = rvblten(1,:,1); end if
            if (present(rthblten)) then; if (allocated(rthblten)) rthblten_loc(1,:,1) = rthblten(1,:,1); end if
            if (present(rqvblten)) then; if (allocated(rqvblten)) rqvblten_loc(1,:,1) = rqvblten(1,:,1); end if
            if (present(rqcblten)) then; if (allocated(rqcblten)) rqcblten_loc(1,:,1) = rqcblten(1,:,1); end if
            if (present(rqiblten)) then; if (allocated(rqiblten)) rqiblten_loc(1,:,1) = rqiblten(1,:,1); end if

            if (present(qc_bl)) then; if (allocated(qc_bl)) qc_bl_loc(1,:,1) = qc_bl(1,:,1); end if
            if (present(qi_bl)) then; if (allocated(qi_bl)) qi_bl_loc(1,:,1) = qi_bl(1,:,1); end if
            if (present(cldfra_bl)) then; if (allocated(cldfra_bl)) cldfra_bl_loc(1,:,1) = cldfra_bl(1,:,1); end if
            if (present(el_pbl)) then; if (allocated(el_pbl)) el_pbl_loc(1,:,1) = el_pbl(1,:,1); end if
            if (present(qke)) then; if (allocated(qke)) qke_loc(1,:,1) = qke(1,:,1); end if
            if (present(tsq)) then; if (allocated(tsq)) tsq_loc(1,:,1) = tsq(1,:,1); end if
            if (present(qsq)) then; if (allocated(qsq)) qsq_loc(1,:,1) = qsq(1,:,1); end if
            if (present(cov)) then; if (allocated(cov)) cov_loc(1,:,1) = cov(1,:,1); end if
            if (present(sh)) then; if (allocated(sh)) sh3d_loc(1,:,1) = sh(1,:,1); end if
            if (present(sm)) then; if (allocated(sm)) sm3d_loc(1,:,1) = sm(1,:,1); end if
            if (present(qke_adv)) then; if (allocated(qke_adv)) qke_adv_loc(1,:,1) = qke_adv(1,:,1); end if
            if (present(pblh)) then; if (allocated(pblh)) pblh_loc(1,1) = pblh(1,1); end if
        end if

        ! Loop through each timestep
        do t = t_start, t_end
            print *, "Processing timestep ", t, " of ", t_end

            ! Read 2D for this timestep: (t,1,1)
            status = nf90_inq_varid(ncid, "XLAND", varid)
            status = nf90_get_var(ncid, varid, xland, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "PSFC", varid)
            status = nf90_get_var(ncid, varid, ps, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "TSK", varid)
            status = nf90_get_var(ncid, varid, ts, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "QSFC", varid)
            status = nf90_get_var(ncid, varid, qsfc, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "UST", varid)
            status = nf90_get_var(ncid, varid, ust, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "FLHC", varid)
            status = nf90_get_var(ncid, varid, ch, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "HFX", varid)
            status = nf90_get_var(ncid, varid, hfx, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "QFX", varid)
            status = nf90_get_var(ncid, varid, qfx, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "WSPD", varid)
            status = nf90_get_var(ncid, varid, wspd, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "ZNT", varid)
            status = nf90_get_var(ncid, varid, znt, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "UOCE", varid)
            status = nf90_get_var(ncid, varid, uoce, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "VOCE", varid)
            status = nf90_get_var(ncid, varid, voce, &
                                  start=[1,1,t], count=[1,1,1])

            dx2d = dxc
            
            ! Read 3D for this timestep:(t,nz,1,1)
            status = nf90_inq_varid(ncid, "dz", varid)
            status = nf90_get_var(ncid, varid, dz_loc(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            ! status = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
            ! print *, "dz has", ndims, "dimensions"

            ! do i = 1, ndims
            !     status = nf90_inquire_dimension(ncid, dimids(i), name=dimname, len=dimlen(i))
            !     print *, "  Dimension", i, ":", trim(dimname), " = ", dimlen(i)
            ! end do

            if (t==2 .and. .not. is_warm_start) then
              initflag = 1
              status = nf90_inq_varid(ncid, "U_mass", varid)
              status = nf90_get_var(ncid, varid, u_loc(1,:,1), &
                                    start=[1,t], count=[nz,1])

              status = nf90_inq_varid(ncid, "V_mass", varid)
              status = nf90_get_var(ncid, varid, v_loc(1,:,1), &
                                    start=[1,t], count=[nz,1])

              status = nf90_inq_varid(ncid, "TH", varid)
              status = nf90_get_var(ncid, varid, th_loc(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])

              status = nf90_inq_varid(ncid, "QVAPOR", varid)
              status = nf90_get_var(ncid, varid, qv_loc(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])

              status = nf90_inq_varid(ncid, "QCLOUD", varid)
              status = nf90_get_var(ncid, varid, qc_loc(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])

              status = nf90_inq_varid(ncid, "QICE", varid)
              status = nf90_get_var(ncid, varid, qi_loc(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])
 
            !else if (t == t_start .and. is_warm_start) then
            else
              initflag = 0
              !print *, 'RUBLTEN(1,:,1)*delt', RUBLTEN_loc(1,:,1)*delt
              u_loc(1,:,1)=u_loc(1,:,1)+RUBLTEN_loc(1,:,1)*delt
              v_loc(1,:,1)=v_loc(1,:,1)+RVBLTEN_loc(1,:,1)*delt
              th_loc(1,:,1)=th_loc(1,:,1)+RTHBLTEN_loc(1,:,1)*delt
              qc_loc(1,:,1)=qc_loc(1,:,1)+RQCBLTEN_loc(1,:,1)*delt
              qv_loc(1,:,1)=qv_loc(1,:,1)+RQVBLTEN_loc(1,:,1)*delt
              qi_loc(1,:,1)=qi_loc(1,:,1)+RQIBLTEN_loc(1,:,1)*delt          
            end if

            status = nf90_inq_varid(ncid, "W_mass", varid)
            status = nf90_get_var(ncid, varid, w_loc(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])


            status = nf90_inq_varid(ncid, "TK", varid)
            status = nf90_get_var(ncid, varid, t3d_loc(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "PTOTAL", varid)
            status = nf90_get_var(ncid, varid, p_loc(1,:,1), &
                 start=[1,1,1,t], count=[1,1,nz,1])

            !calculate p on interfaces:
            pint(:,1,:) = ps
            do k=kts+1,nz
               pint(:,k,:) = p_loc(:,k,:) - p_loc(:,k+1,:)
            enddo
            pint(:,kme,:) = pint(:,kme-1,:)-1.

            status = nf90_inq_varid(ncid, "EXNER", varid)
            status = nf90_get_var(ncid, varid, exner_loc(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "RHO", varid)
            status = nf90_get_var(ncid, varid, rho_loc(1,:,1), &
                                  start=[1,1,1,t+1], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "EXCH_H_mass", varid)
            status = nf90_get_var(ncid, varid, exch_h_loc(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])
            status = nf90_inq_varid(ncid, "EXCH_M_mass", varid)
            status = nf90_get_var(ncid, varid, exch_m_loc(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            ! assign pseudo smoke/chem data
            frp_mean         = 1.0
            emis_ant_no      = 5.0e-9
            vd3d(1, 1, 1)    = 0.01
            vd3d(1, 2:nz, 1) = 1.0e-4
            chem3d           = 1.0e-9
            settle3d         = 0.001


           print *, 'before mynnedmf_driver'

            call mynnedmf_driver    &
                 (ids=ids              , ide=ide             , jds=jds             , jde=jde            , &
                  kds=kds              , kde=kde             , ims=ims             , ime=ime            , &
                  jms=jms              , jme=jme             , kms=kms             , kme=kme            , &
                  its=its              , ite=ite             , jts=jts             , jte=jte            , &
                  kts=kts              , kte=kte             , f_qc=flag_qc        , f_qi=flag_qi       , &
                  f_qs=flag_qs         , f_qnc=flag_qnc      , f_qni=flag_qni      , f_qoz=flag_oz      , &
                  f_qnifa=flag_qnifa   , f_qnwfa=flag_qnwfa  , f_qnbca=flag_qnbca  , initflag=initflag  , &
                  restart=restart1     , cycling=cycling     , delt=delt           ,                      &
                  dx=dx2d              , xland=xland         , ps=ps               , ts=ts              , &
                  qsfc=qsfc            , ust=ust             , ch=ch               , hfx=hfx            , &
                  qfx=qfx              , wspd=wspd           , znt=znt             ,                      &
                  uoce=uoce            , voce=voce           , dz=dz_loc               , u=u_loc                        , &
                  v=v_loc                  , w=w_loc                 , th=th_loc               , tk=t3d_loc             , &
                  p=p_loc                  , pint=pint               , exner=exner_loc         , rho=rho_loc            , &
                  qv=qv_loc                ,                                                                              &
                  qc=qc_loc                , qi=qi_loc               , qs=qs_loc               , qnc=qnc_loc            , &
                  qni=qni_loc              , qnifa=qnifa_loc         , qnwfa=qnwfa_loc         , qnbca=qnbca_loc        , &
!                  qoz=qoz              ,                                                                  &
                  rthraten=rthraten_loc    , pblh=pblh_loc           , kpbl=kpbl           , maxwidth_dd=maxwidth_dd    , &
                  cldfra_bl=cldfra_bl_loc  , qc_bl=qc_bl_loc         , qi_bl=qi_bl_loc          , maxwidth=maxwidth     , &
                  maxmf=maxmf          , ztop_plume=ztop_plume, excess_h=excess_h   , excess_q=excess_q ,                 &
                  maxmf_dd=maxmf_dd    , maxtkeprod=maxtkeprod, cldtop_cooling=cldtop_cooling, ent_eff=ent_eff,           &
                  lwp=lwp              , iwp=iwp              , swp=swp,                                                  &
                  wspd10=wspd10        , wspd80=wspd80        , wspd160=wspd160     , cldceil=cldceil,                    &
                  maxcldfra=maxcldfra  , maxcldfra_pbl=maxcldfra_pbl ,                                                    &
                  qke=qke_loc              , qke_adv=qke_adv_loc     ,                                                    &
                  tsq=tsq_loc              , qsq=qsq_loc             , cov=cov_loc             ,                          &
                  el_pbl=el_pbl_loc        , rublten=rublten_loc     , rvblten=rvblten_loc     , rthblten=rthblten_loc  , &
                  rqvblten=rqvblten_loc    , rqcblten=rqcblten_loc   , rqiblten=rqiblten_loc   , rqsblten=rqsblten_loc  , &
                  rqncblten=rqncblten_loc  , rqniblten=rqniblten_loc , rqnifablten=rqnifablten_loc,rqnwfablten=rqnwfablten_loc, &
                  rqnbcablten=rqnbcablten_loc,                                                                            &
!                  rqozblten=rqozblten  ,                                                                 &
                  edmf_a=edmf_a_loc        , edmf_w=edmf_w_loc       ,                                                    &
                  edmf_qt=edmf_qt_loc      , edmf_thl=edmf_thl_loc   , edmf_ent=edmf_ent_loc   , edmf_qc=edmf_qc_loc    , &
                  sub_thl=sub_thl3d_loc    , sub_sqv=sub_sqv3d_loc   , det_thl=det_thl3d_loc   , det_sqv=det_sqv3d_loc  , &
                  exch_h=exch_h_loc        , exch_m=exch_m_loc       , dqke=dqke_loc           , qwt=qwt_loc            , &
                  qshear=qshear_loc        , qbuoy=qbuoy_loc         , qdiss=qdiss_loc         , sh3d=sh3d_loc          , &
                  sm3d=sm3d_loc            , spp_pbl=spp_pbl     , pattern_spp=pattern_spp_pbl,           &
                  bl_mynn_tkeadvect    = bl_mynn_tkeadvect   , tke_budget          = tke_budget         , &
                  bl_mynn_cloudpdf     = bl_mynn_cloudpdf    , bl_mynn_mixlength   = bl_mynn_mixlength  , &
                  bl_mynn_closure      = bl_mynn_closure     , bl_mynn_edmf        = bl_mynn_edmf       , &
                  bl_mynn_edmf_mom     = bl_mynn_edmf_mom    , bl_mynn_edmf_tke    = bl_mynn_edmf_tke   , &
                  bl_mynn_output       = bl_mynn_output      , bl_mynn_mixscalars  = bl_mynn_mixscalars , &
                  bl_mynn_mixaerosols  = bl_mynn_mixaerosols , bl_mynn_mixnumcon   = bl_mynn_mixnumcon  , &
                  bl_mynn_cloudmix     = bl_mynn_cloudmix    , bl_mynn_mixqt       = bl_mynn_mixqt      , &
                  bl_mynn_edmf_dd      = bl_mynn_edmf_dd     , bl_mynn_ess         = bl_mynn_ess        , &
                  bl_mynn_diags        = bl_mynn_diags       ,                                            &
!#if(WRF_CHEM == 1)
                  mix_chem=mix_chem     , chem3d=chem3d         , vd3d=vd3d     , nchem=nchem           , &
                  ndvel=ndvel           ,                                                                 &
                  settle3d=settle3d     ,                                                                 &
                  frp_mean=frp_mean    , emis_ant_no=emis_ant_no       , enh_mix=enh_mix                , &
!#endif
                  errmsg=errmsg        , errflg=errflg                                                    &
                  )
             
            if (case /= 'clr' .and. bl_mynn_diags >= 2) then
                print '(A, I4, 5(A, F8.3))', 't =', t, ' | U = ', u_loc(1, 1, 1), ' | LWP = ',            &
                        lwp(1, 1), ' | Ceiling = ', cldceil(1,1), ' | Maxcldfra = ', maxcldfra(1,1),      &
                        ' | Maxcldfra_pbl = ', maxcldfra_pbl(1,1)
            else
                print '(A, I4, 1(A, F8.3))', 't =', t, ' | U = ', u_loc(1, 1, 1)
            endif
            
        enddo
        
        ! Close file and deallocate
        status = nf90_close(ncid)

        print *, "Finished processing all timesteps ", t_start, " to ", t_end

        ! Hand updated arrays back to caller if arguments were supplied
        if (present(u)) then
            if (.not. allocated(u)) allocate(u(ims:ime, kms:kme, jms:jme))
            u = u_loc
        end if
        if (present(v)) then
            if (.not. allocated(v)) allocate(v(ims:ime, kms:kme, jms:jme))
            v = v_loc
        end if
        if (present(th)) then
            if (.not. allocated(th)) allocate(th(ims:ime, kms:kme, jms:jme))
            th = th_loc
        end if
        if (present(qv)) then
            if (.not. allocated(qv)) allocate(qv(ims:ime, kms:kme, jms:jme))
            qv = qv_loc
        end if
        if (present(qc)) then
            if (.not. allocated(qc)) allocate(qc(ims:ime, kms:kme, jms:jme))
            qc = qc_loc
        end if
        if (present(qi)) then
            if (.not. allocated(qi)) allocate(qi(ims:ime, kms:kme, jms:jme))
            qi = qi_loc
        end if
        if (present(pblh)) then
            if (.not. allocated(pblh)) allocate(pblh(ims:ime, jms:jme))
            pblh = pblh_loc
        end if
        if (present(rublten)) then
            if (.not. allocated(rublten)) allocate(rublten(ims:ime, kms:kme, jms:jme))
            rublten = rublten_loc
        end if
        if (present(rvblten)) then
            if (.not. allocated(rvblten)) allocate(rvblten(ims:ime, kms:kme, jms:jme))
            rvblten = rvblten_loc
        end if
        if (present(rthblten)) then
            if (.not. allocated(rthblten)) allocate(rthblten(ims:ime, kms:kme, jms:jme))
            rthblten = rthblten_loc
        end if
        if (present(rqvblten)) then
            if (.not. allocated(rqvblten)) allocate(rqvblten(ims:ime, kms:kme, jms:jme))
            rqvblten = rqvblten_loc
        end if
        if (present(rqcblten)) then
            if (.not. allocated(rqcblten)) allocate(rqcblten(ims:ime, kms:kme, jms:jme))
            rqcblten = rqcblten_loc
        end if
        if (present(rqiblten)) then
            if (.not. allocated(rqiblten)) allocate(rqiblten(ims:ime, kms:kme, jms:jme))
            rqiblten = rqiblten_loc
        end if
        if (present(qc_bl)) then
            if (.not. allocated(qc_bl)) allocate(qc_bl(ims:ime, kms:kme, jms:jme))
            qc_bl = qc_bl_loc
        end if
        if (present(qi_bl)) then
            if (.not. allocated(qi_bl)) allocate(qi_bl(ims:ime, kms:kme, jms:jme))
            qi_bl = qi_bl_loc
        end if
        if (present(cldfra_bl)) then
            if (.not. allocated(cldfra_bl)) allocate(cldfra_bl(ims:ime, kms:kme, jms:jme))
            cldfra_bl = cldfra_bl_loc
        end if
        if (present(el_pbl)) then
            if (.not. allocated(el_pbl)) allocate(el_pbl(ims:ime, kms:kme, jms:jme))
            el_pbl = el_pbl_loc
        end if
        if (present(qke)) then
            if (.not. allocated(qke)) allocate(qke(ims:ime, kms:kme, jms:jme))
            qke = qke_loc
        end if
        if (present(qsq)) then
            if (.not. allocated(qsq)) allocate(qsq(ims:ime, kms:kme, jms:jme))
            qsq = qsq_loc
        end if
        if (present(tsq)) then
            if (.not. allocated(tsq)) allocate(tsq(ims:ime, kms:kme, jms:jme))
            tsq = tsq_loc
        end if
        if (present(cov)) then
            if (.not. allocated(cov)) allocate(cov(ims:ime, kms:kme, jms:jme))
            cov = cov_loc
        end if
        if (present(sh)) then
            if (.not. allocated(sh)) allocate(sh(ims:ime, kms:kme, jms:jme))
            sh = sh3d_loc
        end if
        if (present(sm)) then
            if (.not. allocated(sm)) allocate(sm(ims:ime, kms:kme, jms:jme))
            sm = sm3d_loc
        end if
        if (present(qke_adv)) then
            if (.not. allocated(qke_adv)) allocate(qke_adv(ims:ime, kms:kme, jms:jme))
            qke_adv = qke_adv_loc
        end if

        deallocate(time)
        
        ! deallocate 2D arrays
        deallocate(xland,ps,ts,qsfc,ust,ch,hfx,qfx,wspd,znt,uoce,voce,                &
             kpbl,maxmf,maxwidth,pblh_loc,ztop_plume,excess_h,excess_q,               &
             maxwidth_dd,maxmf_dd,maxtkeprod,cldtop_cooling,ent_eff)
        
        ! deallocate 3D arrays
        deallocate(u_loc, v_loc, w_loc, th_loc, t3d_loc, p_loc, exner_loc, rho_loc,   &
                    qv_loc, qc_loc, qi_loc, dz_loc, exch_h_loc, exch_m_loc, pint)
        deallocate(cov_loc, det_thl3d_loc, det_sqv3d_loc, dqke_loc, edmf_a_loc,       &
                    edmf_ent_loc, edmf_qc_loc, edmf_qt_loc, edmf_thl_loc, edmf_w_loc, &
                    qnc_loc, qni_loc, qnwfa_loc, qnifa_loc, qs_loc, qshear_loc,       &
                    qBUOY_loc, qDISS_loc) 
        deallocate(qwt_loc, qke_loc, qke_adv_loc, tsq_loc, qsq_loc, el_pbl_loc,       &
                    sh3d_loc, sm3d_loc, qc_bl_loc, qi_bl_loc,                         &
                   cldfra_bl_loc, sub_thl3d_loc, sub_sqv3d_loc,                       & 
                   RQVBLTEN_loc, RQCBLTEN_loc, RQIBLTEN_loc)
                  ! , qnbca_loc, &
                  ! rqsblten_loc, rqncblten_loc, rqniblten_loc, rqnifablten_loc, &
                  ! rqnwfablten_loc, rqnbcablten_loc)
   end subroutine mynnedmf_test


   subroutine restart_test(case,bl_mynn_closure,bl_mynn_cloudpdf,bl_mynn_mixlength,   &
                                bl_mynn_edmf,bl_mynn_edmf_dd,bl_mynn_edmf_mom,        &
                                bl_mynn_edmf_tke,bl_mynn_cloudmix,                    &
                                bl_mynn_mixqt, bl_mynn_mixscalars,                    &
                                bl_mynn_mixaerosols,bl_mynn_mixnumcon,                &
                                bl_mynn_ess,tke_budget,bl_mynn_diags,                 &
                                n_restart_in,tol_in)

        implicit none
        character(len=*),intent(in)  :: case
        integer, intent(in)          :: bl_mynn_cloudpdf, bl_mynn_mixlength, bl_mynn_edmf,   &
                                        bl_mynn_edmf_dd, bl_mynn_edmf_mom, bl_mynn_edmf_tke, &
                                        bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_mixscalars, &
                                        bl_mynn_mixaerosols, bl_mynn_mixnumcon, bl_mynn_ess, & 
                                        tke_budget, bl_mynn_diags
        real, intent(in)             :: bl_mynn_closure
        integer, intent(in), optional :: n_restart_in
        real,    intent(in), optional :: tol_in

        integer :: n_restart
        real    :: tol
        logical :: test_pass
        real    :: max_du, max_dth, max_dqv, max_dqke, max_dpblh

        ! baseline cold start
        real, dimension(:,:,:), allocatable :: u_base,v_base,th_base,qv_base,qc_base,qi_base,qke_base
        real, dimension(:,:), allocatable :: pblh_base

        ! before restart
        real, dimension(:,:,:), allocatable :: u1,v1,th1,qv1,qc1,qi1
        real, dimension(:,:,:), allocatable :: ru1,rv1,rth1,rqv1,rqc1,rqi1
        real, dimension(:,:,:), allocatable :: qc_bl1,qi_bl1,cldfra_bl1
        real, dimension(:,:,:), allocatable :: el_pbl1,qke1,qsq1,tsq1,cov1,sh1,sm1,qke_adv1
        real, dimension(:,:), allocatable :: pblh1
        ! after restart
        real, dimension(:,:,:), allocatable :: u2,v2,th2,qv2,qc2,qi2,qke2
        real, dimension(:,:), allocatable :: pblh2

        n_restart = 10
        if (present(n_restart_in)) n_restart = n_restart_in
        tol = 1.0e-6
        if (present(tol_in)) tol = tol_in

        print*, '=== MYNN-EDMF restart test, case: ', trim(case), ', n_restart = ', n_restart, ' ==='

        !--- before restart: t=2..n_restart ---
        call mynnedmf_test(case=case,bl_mynn_closure=bl_mynn_closure,bl_mynn_cloudpdf=bl_mynn_cloudpdf,                &
                bl_mynn_mixlength=bl_mynn_mixlength,bl_mynn_edmf=bl_mynn_edmf,bl_mynn_edmf_dd=bl_mynn_edmf_dd,         &
                bl_mynn_edmf_mom=bl_mynn_edmf_mom,bl_mynn_edmf_tke=bl_mynn_edmf_tke,bl_mynn_cloudmix=bl_mynn_cloudmix, &
                bl_mynn_mixqt=bl_mynn_mixqt,bl_mynn_mixscalars=bl_mynn_mixscalars,                                     &
                bl_mynn_mixaerosols=bl_mynn_mixaerosols,bl_mynn_mixnumcon=bl_mynn_mixnumcon,                           &
                bl_mynn_ess=bl_mynn_ess,tke_budget=tke_budget,bl_mynn_diags=bl_mynn_diags,                             &
                t_end_in=n_restart,                                                                                    &
                u=u1,v=v1,th=th1,qv=qv1,qc=qc1,qi=qi1,                                                                 &
                rublten=ru1,rvblten=rv1,rthblten=rth1,                                                                 &
                rqvblten=rqv1,rqcblten=rqc1,rqiblten=rqi1,                                                             &
                qc_bl=qc_bl1,qi_bl=qi_bl1,cldfra_bl=cldfra_bl1,                                                        &
                el_pbl=el_pbl1,qke=qke1,tsq=tsq1,qsq=qsq1,cov=cov1,                                                    &
                sh=sh1,sm=sm1,qke_adv=qke_adv1,pblh=pblh1)

        ! --- after restart, t=n_restart+1..nt---
        call mynnedmf_test(case=case,bl_mynn_closure=bl_mynn_closure,bl_mynn_cloudpdf=bl_mynn_cloudpdf,                &
                bl_mynn_mixlength=bl_mynn_mixlength,bl_mynn_edmf=bl_mynn_edmf,bl_mynn_edmf_dd=bl_mynn_edmf_dd,         &
                bl_mynn_edmf_mom=bl_mynn_edmf_mom,bl_mynn_edmf_tke=bl_mynn_edmf_tke,bl_mynn_cloudmix=bl_mynn_cloudmix, &
                bl_mynn_mixqt=bl_mynn_mixqt,bl_mynn_mixscalars=bl_mynn_mixscalars,                                     &
                bl_mynn_mixaerosols=bl_mynn_mixaerosols,bl_mynn_mixnumcon=bl_mynn_mixnumcon,                           &
                bl_mynn_ess=bl_mynn_ess,tke_budget=tke_budget,bl_mynn_diags=bl_mynn_diags,                             &
                t_start_in=n_restart+1, restart_in=.true.,                                                             &
                u=u1,v=v1,th=th1,qv=qv1,qc=qc1,qi=qi1,                                                                 &
                rublten=ru1,rvblten=rv1,rthblten=rth1,                                                                 &
                rqvblten=rqv1,rqcblten=rqc1,rqiblten=rqi1,                                                             &
                qc_bl=qc_bl1,qi_bl=qi_bl1,cldfra_bl=cldfra_bl1,                                                        &
                el_pbl=el_pbl1,qke=qke1,qsq=qsq1,tsq=tsq1,                                                             &
                cov=cov1,sh=sh1,sm=sm1,qke_adv=qke_adv1,                                                               &
                pblh=pblh1)

        !--- baseline: continuous run, t=2..nt ---
        call mynnedmf_test(case=case,bl_mynn_closure=bl_mynn_closure,bl_mynn_cloudpdf=bl_mynn_cloudpdf,                &
                bl_mynn_mixlength=bl_mynn_mixlength,bl_mynn_edmf=bl_mynn_edmf,bl_mynn_edmf_dd=bl_mynn_edmf_dd,         &
                bl_mynn_edmf_mom=bl_mynn_edmf_mom,bl_mynn_edmf_tke=bl_mynn_edmf_tke,bl_mynn_cloudmix=bl_mynn_cloudmix, &
                bl_mynn_mixqt=bl_mynn_mixqt,bl_mynn_mixscalars=bl_mynn_mixscalars,                                     &
                bl_mynn_mixaerosols=bl_mynn_mixaerosols,bl_mynn_mixnumcon=bl_mynn_mixnumcon,                           &
                bl_mynn_ess=bl_mynn_ess,tke_budget=tke_budget,bl_mynn_diags=bl_mynn_diags,                             &
                restart_in=.false.,                                                                                    &
                u=u_base,v=v_base,th=th_base,qv=qv_base,qc=qc_base,qi=qi_base,                                         &
                qke=qke_base,pblh=pblh_base)

        ! --- compare baseline vs restarted final state ---
        max_du    = maxval(abs(u_base    - u1))
        max_dth   = maxval(abs(th_base   - th1))
        max_dqv   = maxval(abs(qv_base   - qv1))
        max_dqke  = maxval(abs(qke_base  - qke1))
        max_dpblh = maxval(abs(pblh_base - pblh1))

        test_pass = (max_du < tol) .and. (max_dth < tol) .and. (max_dqv < tol) .and. &
                        (max_dqke < tol) .and. (max_dpblh < tol)

        print *, '--- MYNN-EDMF restart test results, case: ', trim(case), ' ---'
        print *, '  max|du|   =', max_du
        print *, '  max|dth|  =', max_dth
        print *, '  max|dqv|  =', max_dqv
        print *, '  max|dqke| =', max_dqke
        print *, '  max|dpblh|=', max_dpblh
        if (test_pass) then
                print *, '  RESULT: Restart Test PASS (tol=', tol, ')'
        else
                print *, '  RESULT: Restart Test FAIL (tol=', tol, ')'
                stop 1
        end if

   end subroutine restart_test

end module module_bl_mynnedmf_tests

