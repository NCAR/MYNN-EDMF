!> \file module_bl_mynnedmf_driver.F90
!!  This serves as the unified interface between the WRF/MPAS PBL drivers and the MYNN-EDMF
!!  turbulence scheme (in module_bl_mynnedmf.F90).

!>\ingroup gsd_mynn_edmf
!> The following references best describe the code within
!!    Olson et al. (2026, NOAA Technical Memorandum)
!!    Nakanishi and Niino (2009) \cite NAKANISHI_2009
!=================================================================================================================
 module module_bl_mynnedmf_driver

 use module_bl_mynnedmf_diags, only: mynnedmf_diags
 use module_bl_mynnedmf_common,only: kind_phys,xlvcp,xlscp
 use module_bl_mynnedmf,only: mynnedmf

 implicit none
 real(kind_phys),parameter::zero=0.0
 real(kind_phys),parameter::one =1.0
 
 private
 public:: mynnedmf_driver
 public:: mynnedmf_init
 public:: mynnedmf_finalize

 contains

!=================================================================================================================
!> \section arg_table_mynnedmf_init Argument Table                                                                                                                           
!! \htmlinclude mynnedmf_wrapper_init.html                                                                                                                                   
!!                                                                                                                                                                           
 subroutine mynnedmf_init (                        &
   &  rublten,rvblten,rthblten,rqvblten,rqcblten,  &
   &  rqiblten,qke,                                &
   &  restart,allowed_to_read,                     &
   &  p_qc,p_qi,param_first_scalar,                &
   &  ids,ide,jds,jde,kds,kde,                     &
   &  ims,ime,jms,jme,kms,kme,                     &
   &  its,ite,jts,jte,kts,kte                      )

   use module_bl_mynnedmf_common, only: kind_phys,zero

   implicit none

   logical,intent(in) :: allowed_to_read,restart

   integer,intent(in) :: ids,ide,jds,jde,kds,kde,  &                                                                                                                         
        &                ims,ime,jms,jme,kms,kme,  &
        &                its,ite,jts,jte,kts,kte

   real(kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout) :: &                                                                                                                
        &rublten,rvblten,rthblten,rqvblten,                 &
        &rqcblten,rqiblten,qke

   integer,  intent(in) :: p_qc,p_qi,param_first_scalar

   integer :: i,j,k,itf,jtf,ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   if (.not.restart) then
      do j=jts,jtf
      do k=kts,ktf
      do i=its,itf
         rublten(i,k,j)=zero
         rvblten(i,k,j)=zero
         rthblten(i,k,j)=zero
         rqvblten(i,k,j)=zero
         if( p_qc >= param_first_scalar ) rqcblten(i,k,j)=zero
         if( p_qi >= param_first_scalar ) rqiblten(i,k,j)=zero
      enddo
      enddo
      enddo
   endif

 end subroutine mynnedmf_init
!=================================================================================================================
 subroutine mynnedmf_finalize ()
 end subroutine mynnedmf_finalize
!================================================================================================================= 
 SUBROUTINE mynnedmf_driver(           &
                  ids               , ide               , jds                , jde                , &
                  kds               , kde               , ims                , ime                , &
                  jms               , jme               , kms                , kme                , &
                  its               , ite               , jts                , jte                , &
                  kts               , kte               , f_qc               , f_qi               , &
                  f_qs              , f_qoz             , f_qnc              , f_qni              , &
                  f_qnifa           , f_qnwfa           , f_qnbca            , initflag           , &
                  restart           , cycling           , delt               ,                      &
                  dx                , xland             , ps                 , ts                 , &
                  qsfc              , ust               , ch                 , hfx                , &
                  qfx               , wspd              , znt                ,                      &
                  uoce              , voce              ,                                           &
                  !3d input
                  dz                , u                 , v                  , w                  , &
                  th                , tk                , p                  , exner              , &
                  rho               , qv                , qc                 , qi                 , &
                  qs                , qnc               , qni                , qnifa              , &
                  qnwfa             , qnbca             , qoz                , rthraten           , &
                  pint              ,                                                               &
                  !3d output
                  cldfra_bl         , qc_bl             , qi_bl              ,                      &
                  qke               , qke_adv           ,                                           &
                  tsq               , qsq               , cov                , el_pbl             , &
                  !output tendencies
                  rublten           , rvblten           , rthblten           ,                      &
                  rqvblten          , rqcblten          , rqiblten           , rqsblten           , &
                  rqncblten         , rqniblten         , rqnifablten        , rqnwfablten        , &
                  rqnbcablten       , rqozblten         ,                                           &
                  !2d output
                  pblh              , kpbl              , maxwidth           ,                      &
                  maxmf             , ztop_plume        , excess_h           , excess_q           , &
                  maxwidth_dd       , maxmf_dd          , maxtkeprod         , cldtop_cooling     , &
                  ent_eff           ,                                                               &
                  !optional 2d diagnostic output
                  lwp               , iwp               , swp                , wspd10             , &
                  wspd80            , wspd160           , cldceil            , maxcldfra          , &
                  maxcldfra_pbl     ,                                                               &
                  !optional 3d output
                  edmf_a            , edmf_w            ,                                           &
                  edmf_qt           , edmf_thl          , edmf_ent           , edmf_qc            , &
                  sub_thl           , sub_sqv           , det_thl            , det_sqv            , &
                  exch_h            , exch_m            , dqke               , qwt                , &
                  qshear            , qbuoy             , qdiss              , sh3d               , &
                  sm3d              ,                                                               &
                  !configuration options (+spp array)
                  spp_pbl           , pattern_spp       ,                                           &
                  bl_mynn_tkeadvect , tke_budget        , bl_mynn_cloudpdf   , bl_mynn_mixlength  , &
                  bl_mynn_closure   , bl_mynn_edmf      , bl_mynn_edmf_mom   , bl_mynn_edmf_tke   , &
                  bl_mynn_output    , bl_mynn_mixscalars, bl_mynn_mixaerosols, bl_mynn_mixnumcon  , &
                  bl_mynn_cloudmix  , bl_mynn_mixqt     , bl_mynn_edmf_dd    , bl_mynn_ess        , &
                  bl_mynn_diags     ,                                                               &
                  !smoke/dust
                  mix_chem          , nchem             , ndvel              , enh_mix            , &
                  chem3d            , settle3d          , vd3d               ,                      &
                  frp_mean          , emis_ant_no       ,                                           &
                  !ccpp error handling
                  errmsg            , errflg                                                        &
               )

!=================================================================================================================
 use module_bl_mynnedmf_common,only: kind_phys,zero,one
   
!--- input arguments:
 logical,intent(in):: &
    f_qc,               &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi,               &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs,               &! if true,the physics package includes the snow mixing ratio.
    f_qoz,              &! if true,the physics package includes the ozone mixing ratio.
    f_qnc,              &! if true,the physics package includes the cloud liquid water number concentration.
    f_qni,              &! if true,the physics package includes the cloud ice number concentration.
    f_qnifa,            &! if true,the physics package includes the "ice-friendly" aerosol number concentration.
    f_qnwfa,            &! if true,the physics package includes the "water-friendly" aerosol number concentration.
    f_qnbca               ! if true,the physics package includes the number concentration of black carbon.

 logical,intent(in):: &
    bl_mynn_tkeadvect,  &
    restart,            &!
    cycling

 integer,intent(in):: &
    ids,ide,jds,jde,kds,kde, &
    ims,ime,jms,jme,kms,kme, &
    its,ite,jts,jte,kts,kte

 integer,intent(in):: &
    bl_mynn_cloudpdf,   &!
    bl_mynn_mixlength,  &!
    bl_mynn_edmf,       &!
    bl_mynn_edmf_dd,    &!
    bl_mynn_edmf_mom,   &!
    bl_mynn_edmf_tke,   &!
    bl_mynn_output,     &!
    bl_mynn_mixscalars, &!
    bl_mynn_mixaerosols,&!
    bl_mynn_mixnumcon,  &!
    bl_mynn_cloudmix,   &!
    bl_mynn_mixqt,      &!
    bl_mynn_ess,        &!
    tke_budget,         &!
    bl_mynn_diags
 
 integer,intent(in):: &
    initflag,           &!
    spp_pbl              !

 real(kind_phys),intent(in):: &
    bl_mynn_closure

 real(kind_phys),intent(in):: &
    delt                 !

 real(kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    dx,                 &!
    xland,              &!
    ps,                 &!
    ts,                 &!
    qsfc,               &!
    ust,                &!
    ch,                 &!
    hfx,                &!
    qfx,                &!
    wspd,               &!
    uoce,               &!
    voce,               &!
    znt                  !

 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz,          &!
    u,           &!
    w,           &!
    v,           &!
    th,          &!
    tk,          &!
    p,           &!
    pint,        &!
    exner,       &!
    rho,         &!
    qv,          &!
    rthraten      !

 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    qc,          &!
    qi,          &!
    qs,          &!
    qoz,         &!
    qnc,         &!
    qni,         &!
    qnifa,       &!
    qnwfa,       &!
    qnbca

 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    pattern_spp   !


!--- inout arguments:
 integer,intent(inout),dimension(ims:ime,jms:jme):: &
    kpbl

 real(kind_phys),intent(inout),dimension(ims:ime,jms:jme):: &
    pblh          !

 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    cldfra_bl,   &!
    qc_bl,       &!
    qi_bl         !

 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    el_pbl,      &!
    qke,         &!
    qke_adv,     &!
    cov,         &!
    qsq,         &!
    tsq,         &!
    sh3d,        &!
    sm3d

 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    rublten,     &!
    rvblten,     &!
    rthblten

 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    rqvblten,    &!
    rqcblten,    &!
    rqiblten,    &!
    rqsblten,    &!
    rqozblten,   &!
    rqncblten,   &!
    rqniblten,   &!
    rqnifablten, &!
    rqnwfablten, &!
    rqnbcablten   !

 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    edmf_a,      &!
    edmf_w,      &!
    edmf_qt,     &!
    edmf_thl,    &!
    edmf_ent,    &!
    edmf_qc,     &!
    sub_thl,     &!
    sub_sqv,     &!
    det_thl,     &!
    det_sqv       !

 real(kind_phys),intent(inout),dimension(ims:ime,jms:jme),optional:: &
    lwp,      &!
    iwp,      &!
    swp,      &!
    cldceil
 real(kind_phys),intent(out),dimension(ims:ime,jms:jme),optional::   &
    wspd10,      &!
    wspd80,      &!
    wspd160,     &!
    maxcldfra,   &!
    maxcldfra_pbl

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg        ! output error message (-).

 integer,intent(out):: &
    errflg        ! output error flag (-).

 real(kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    maxwidth,    &!
    maxmf,       &!
    ztop_plume,  &!
    excess_h,    &!
    excess_q,    &
    maxwidth_dd, &
    maxmf_dd,    &
    maxtkeprod,  &
    cldtop_cooling,&
    ent_eff 

 real(kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    exch_h,      &!
    exch_m        !

 real(kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    dqke,        &!
    qwt,         &!
    qshear,      &!
    qbuoy,       &!
    qdiss         !

!--- input arguments for PBL and free-tropospheric mixing of chemical species:
 logical,intent(in),optional:: mix_chem
 logical::mix_chem1
 integer,intent(in),optional:: nchem,ndvel
 integer::nchem1,ndvel1
 logical,intent(in),optional:: enh_mix
 logical::enh_mix1 
! real(kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: frp_mean,emis_ant_no
! real(kind_phys),intent(in),dimension(ims:ime,jms:jme,ndvel),optional:: vd3d
! real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme,nchem),optional:: chem3d,settle3d
 real(kind_phys),intent(in),dimension(:,:),  optional:: frp_mean,emis_ant_no
 real(kind_phys),intent(in),dimension(:,:,:),optional:: vd3d
 real(kind_phys),intent(inout),dimension(:,:,:,:),optional:: chem3d,settle3d
 !local smoke/chem arrays
 real(kind_phys):: frp1,emisant_no1
 !real(kind_phys),dimension(ndvel):: vd1
 !real(kind_phys),dimension(kts:kte,nchem):: chem1,settle1
 real(kind_phys),allocatable,dimension(:):: vd1
 real(kind_phys),allocatable,dimension(:,:):: chem1,settle1
 
!generic scalar array support
 integer, parameter :: nscalars=1
 real(kind_phys),dimension(kts:kte,nscalars):: scalars

!local
 integer:: i,k,j,nc,nd

 integer:: kpbl1

 real(kind_phys):: max_hfx

 real(kind_phys):: &
    dx1,xland1,ps1,ts1,qsfc1,ust1,ch1,hfx1,qfx1, &
    wspd1,uoce1,voce1,znt1

 real(kind_phys),dimension(kts:kte):: &
    dz1,u1,v1,th1,tk1,p1,exner1,rho1,qv1,rthraten1,delp1,zagl1

 real(kind_phys),dimension(kts:kte+1):: zw1

 real(kind_phys),dimension(kts:kme):: &
    w1

 real(kind_phys),dimension(kts:kte):: &
    qc1,qi1,qs1,qnc1,qni1,qnifa1,qnwfa1,qnbca1,qoz1

 real(kind_phys),dimension(kts:kte):: &
    pattern_spp1

 real(kind_phys):: &
    pblh1, lwp1, iwp1, swp1, wspd101, wspd801, wspd1601, cldceil1,      &
    maxcldfra1 , maxcldfra_pbl1

 real(kind_phys),dimension(kts:kte):: &
    cldfra_bl1,qc_bl1,qi_bl1,el_pbl1,qke1,qke_adv1,cov1,qsq1,tsq1,sh1,sm1

 real(kind_phys),dimension(kts:kte):: &
    rublten1,rvblten1,rthblten1,rqvblten1,rqcblten1,rqiblten1,rqsblten1, &
    rqncblten1,rqniblten1,rqnifablten1,rqnwfablten1,rqnbcablten1,rqozblten1

 real(kind_phys),dimension(kts:kte):: &
    edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1, &
    sub_thl1,sub_sqv1,det_thl1,det_sqv1

 real(kind_phys):: &
    maxwidth1,maxmf1,ztop_plume1,excess_h1,excess_q1, &
    maxwidth_dd1,maxmf_dd1,maxtkeprod1,cldtop_cooling1,ent_eff1

 real(kind_phys),dimension(kts:kte):: &
    kh1,km1

 real(kind_phys), dimension(:), allocatable :: &
    dqke1,qWT1,qSHEAR1,qBUOY1,qDISS1

 real(kind_phys),dimension(kts:kte):: &
    sqv1,sqc1,sqi1,sqs1

!-----------------------------------------------------------------------------------------------------------------
 !for debug printing, set to true
 logical, parameter:: debug = .false.

#if ((EM_CORE == 1) || defined(mpas))
 logical, parameter:: dry_mixing_ratio = .true.
 logical, parameter:: neg_moist_check  = .false.
#else
 logical, parameter:: dry_mixing_ratio = .false.
 logical, parameter:: neg_moist_check  = .true.
#endif

 if (debug) then
    write(0,*)"=============================================="
    write(0,*)"dry_mixing_ratio ", dry_mixing_ratio
    write(0,*)"in mynn-edmf driver..."
    write(0,*)"initflag=",initflag," restart =",restart
 endif

 errmsg = " "
 errflg = 0

 !tke budget (optional arrays)
 if (tke_budget == 1) then
    allocate(dqke1(kts:kte),    source=zero)
    allocate(qwt1(kts:kte),     source=zero)
    allocate(qshear1(kts:kte),  source=zero)
    allocate(qbuoy1(kts:kte),   source=zero)
    allocate(qdiss1(kts:kte),   source=zero)
 endif

 if (present(enh_mix)) then
    enh_mix1 = enh_mix
 else
    enh_mix1 = .false.
 endif

 if(present(mix_chem) .and. present(chem3d).and. present(settle3d) .and. present(vd3d) .and. &
    present(frp_mean) .and. present(emis_ant_no)) then
    mix_chem1 = mix_chem
 else
    mix_chem1 = .false.
 endif

 !---------------------------------------
 !Begin looping in the i- and j-direction
 !---------------------------------------
 do j = jts,jte
 do i = its,ite
     
    !--- input arguments
    dx1    = dx(i,j)
    xland1 = xland(i,j)
    ps1    = ps(i,j)
    ts1    = ts(i,j)
    qsfc1  = qsfc(i,j)
    ust1   = ust(i,j)
    ch1    = ch(i,j)
    wspd1  = wspd(i,j)
    uoce1  = uoce(i,j)
    voce1  = voce(i,j)
    znt1   = znt(i,j)
    !check for unearthly incoming surface fluxes. These threshold are only surpassed
    !when something unphysical is happening (or there is a fire). In this unphysical state,
    !conservation is already questionable and the simulation is heading off the rails. 
    !Try to curb the consequences of this behavior by imposing liberal limits on
    !the incoming fluxes:
    max_hfx = 1200._kind_phys
    if (present(frp_mean)) then
       if (frp_mean(i,j) > 10._kind_phys) max_hfx = 6000._kind_phys
    endif
    hfx1 = hfx(i,j)
    if (hfx1 > max_hfx) then
       !print*,"hfx at i=",i," j=",j,"is unrealistic:",hfx1
       hfx1 = max_hfx
    endif
    if (hfx1 < -600.) then
       !print*,"hfx at i=",i," j=",j,"is unrealistic:",hfx1
       hfx1 = -600._kind_phys
    endif
    qfx1 = qfx(i,j)
    if (qfx1 > 9e-4) then
       !print*,"qfx at i=",i," j=",j,"is unrealistic:",qfx1
       qfx1 = 9e-4_kind_phys
    endif
    if (qfx1 < -3e-4) then
       !print*,"qfx at i=",i," j=",j,"is unrealistic:",qfx1
       qfx1 = -3e-4_kind_phys
    endif

    zw1(kts)        = zero
    zagl1(kts)      = zero
    do k = kts,kte
       dz1(k)       = dz(i,k,j)
       u1(k)        = u(i,k,j)
       v1(k)        = v(i,k,j)
       w1(k)        = w(i,k,j)
       th1(k)       = th(i,k,j)
       tk1(k)       = tk(i,k,j)
       p1(k)        = p(i,k,j)
       exner1(k)    = exner(i,k,j)
       rho1(k)      = rho(i,k,j)
       rthraten1(k) = rthraten(i,k,j)
       delp1(k)     = max(0.01_kind_phys, pint(i,k,j)-pint(i,k+1,j))
       zagl1(k)     = zw1(k) + 0.5_kind_phys * dz1(k)                   ! at mass point AGL
       zw1(k+1)     = zw1(k) + dz1(k)                                   ! at interface
    enddo
    w1(kte+1) = w(i,kte+1,j)

    !--- input arguments for cloud mixing ratios and number concentrations; input argument
    !    for the ozone mixing ratio; input arguments for aerosols from the aerosol-aware
    !    Thompson cloud microphysics:
    do k = kts,kte
       qv1(k)    = max(1e-10_kind_phys, qv(i,k,j))
       qc1(k)    = zero
       qi1(k)    = zero
       qs1(k)    = zero
       qoz1(k)   = zero
       qnc1(k)   = zero
       qni1(k)   = zero
       qnifa1(k) = zero
       qnwfa1(k) = zero
       qnbca1(k) = zero
    enddo
    if(f_qc .and. present(qc)) then
       do k = kts,kte
          qc1(k) = qc(i,k,j)
       enddo
    endif
    if(f_qi .and. present(qi)) then
       do k = kts,kte
          qi1(k) = qi(i,k,j)
       enddo
    endif
    if(f_qs .and. present(qs)) then
       do k = kts,kte
          qs1(k) = qs(i,k,j)
       enddo
    endif
    if(f_qnc .and. present(qnc)) then
       do k = kts,kte
          qnc1(k) = qnc(i,k,j)
       enddo
    endif
    if(f_qni .and. present(qni)) then
       do k = kts,kte
          qni1(k) = qni(i,k,j)
       enddo
    endif
    if(f_qnifa .and. present(qnifa)) then
       do k = kts,kte
          qnifa1(k) = qnifa(i,k,j)
       enddo
    endif
    if(f_qnwfa .and. present(qnwfa)) then
       do k = kts,kte
          qnwfa1(k) = qnwfa(i,k,j)
       enddo
    endif
    if(f_qnbca .and. present(qnbca)) then
       do k = kts,kte
          qnbca1(k) = qnbca(i,k,j)
       enddo
    endif
    if(f_qoz .and. present(qoz)) then
       do k = kts,kte
          qoz1(k) = qoz(i,k,j)
       enddo
    endif

    do k = kts,kte
       rqcblten1(k)    = zero
       rqiblten1(k)    = zero
       rqsblten1(k)    = zero
       rqozblten1(k)   = zero
       rqncblten1(k)   = zero
       rqniblten1(k)   = zero
       rqnifablten1(k) = zero
       rqnwfablten1(k) = zero
       rqnbcablten1(k) = zero
    enddo
    
    if (neg_moist_check .and. present(qi) .and. present(qc)) then
       !find/fix negative mixing ratios
       call moisture_check2(kte        , delt       , delp1      , exner1     , &
                            qv1        , qc1        , qi1        , th1          )
       do k = kts,kte
          tk1(k)    = th1(k)*exner1(k)
       enddo
    endif
    
    !--- conversion from mixing ratios to specific contents:
    if (dry_mixing_ratio) then
       call mynnedmf_pre_run(kte,f_qc,f_qi,f_qs,qv1,qc1,qi1,qs1,sqv1,sqc1, &
                            sqi1,sqs1,errmsg,errflg)
    endif

    !--- initialization of the stochastic forcing in the PBL:
    if(spp_pbl > 0 .and. present(pattern_spp)) then
       do k = kts,kte
          pattern_spp1(k) = pattern_spp(i,k,j)
       enddo
    else
       do k = kts,kte
          pattern_spp1(k) = zero
       enddo
    endif

    !--- inout arguments:
    pblh1 = pblh(i,j)
    kpbl1 = kpbl(i,j)

    !when NOT cold-starting on the first time step, update input
    if (initflag .eq. 0 .or. restart) then
      !update sgs cloud info.
      do k = kts,kte
         cldfra_bl1(k) = cldfra_bl(i,k,j)
         qc_bl1(k)     = qc_bl(i,k,j)
         qi_bl1(k)     = qi_bl(i,k,j)
      enddo

      !turbulence variables
      do k = kts,kte
         el_pbl1(k)  = el_pbl(i,k,j)
         qke1(k)    = qke(i,k,j)
      !   qke_adv1(k) = qke_adv(i,k,j)
         cov1(k)    = cov(i,k,j)
         tsq1(k)    = tsq(i,k,j)   
         qsq1(k)    = qsq(i,k,j)
         sh1(k)     = sh3d(i,k,j)
         sm1(k)     = sm3d(i,k,j)
         if (bl_mynn_tkeadvect) then
            qke_adv1(kts:kte) = qke_adv(i,kts:kte,j)
         else
            qke_adv1(kts:kte) = qke(i,kts:kte,j)
         endif          
      enddo
    endif
    
    !Smoke/dust
    if(mix_chem1) then
       ndvel1    = ndvel
       nchem1    = nchem
       allocate(vd1(ndvel1))
       allocate(chem1(kts:kte,nchem1))
       allocate(settle1(kts:kte,nchem1))
       do nc = 1,nchem
         do k = kts,kte
            chem1(k,nc)   = chem3d(i,k,j,nc)
            settle1(k,nc) = settle3d(i,k,j,nc)
         enddo
       enddo
       do nd = 1,ndvel
          vd1(nd) = vd3d(i,j,nd)
       enddo
       frp1        = frp_mean(i,j)
       emisant_no1 = emis_ant_no(i,j)
    else
       ndvel1      = 1
       nchem1      = 1
       allocate(vd1(ndvel1))
       allocate(chem1(kts:kte,nchem1))
       allocate(settle1(kts:kte,nchem1))
       chem1       = zero
       settle1     = zero
       vd1         = zero
       frp1        = zero
       emisant_no1 = zero
    endif
    scalars     = zero


    if (debug) then
      print*,"In mynnedmf driver, just before the call to mynnedmf"
    endif
    
    call mynnedmf( &
            i               = i             , j           = j             ,                              &
            initflag        = initflag      , restart     = restart       , cycling     = cycling      , &
            delt            = delt          , dz1         = dz1           , dx          = dx1          , &
            znt             = znt1          , u1          = u1            , v1          = v1           , &
            w1              = w1            , th1         = th1           , sqv1        = sqv1         , &
            sqc1            = sqc1          , sqi1        = sqi1          , sqs1        = sqs1         , &
            qnc1            = qnc1          , qni1        = qni1          , qnwfa1      = qnwfa1       , &
            qnifa1          = qnifa1        , qnbca1      = qnbca1        , ozone1      = qoz1         , &
            delp1           = delp1         , zw1         = zw1           , zagl1       = zagl1        , &
            pres1           = p1            , ex1         = exner1        , rho1        = rho1         , &
            tk1             = tk1           , xland       = xland1        , ts          = ts1          , &
            qsfc            = qsfc1         , ps          = ps1           , ust         = ust1         , &
            ch              = ch1           , hfx         = hfx1          , qfx         = qfx1         , &
            wspd            = wspd1         , uoce        = uoce1         , voce        = voce1        , &
            qke1            = qke1          , qke_adv1    = qke_adv1      ,                              &
            tsq1            = tsq1          , qsq1        = qsq1          , cov1        = cov1         , &
            rthraten1       = rthraten1     , du1         = rublten1      , dv1         = rvblten1     , &
            dth1            = rthblten1     , dqv1        = rqvblten1     , dqc1        = rqcblten1    , &
            dqi1            = rqiblten1     , dqs1        = rqsblten1     , dqnc1       = rqncblten1   , &
            dqni1           = rqniblten1    , dqnwfa1     = rqnwfablten1  , dqnifa1     = rqnifablten1 , &
            dqnbca1         = rqnbcablten1  , dozone1     = rqozblten1    , kh1         = kh1          , &
            km1             = km1           , pblh        = pblh1         , kpbl        = kpbl1        , &
            el1             = el_pbl1       , dqke1       = dqke1         , qwt1        = qwt1         , &
            qshear1         = qshear1       , qbuoy1      = qbuoy1        , qdiss1      = qdiss1       , &
            sh1             = sh1           , sm1         = sm1           , qc_bl1      = qc_bl1       , &
            qi_bl1          = qi_bl1        , cldfra_bl1  = cldfra_bl1    ,                              &
            edmf_a1         = edmf_a1       , edmf_w1     = edmf_w1       , edmf_qt1    = edmf_qt1     , &
            edmf_thl1       = edmf_thl1     , edmf_ent1   = edmf_ent1     , edmf_qc1    = edmf_qc1     , &
            sub_thl1        = sub_thl1      , sub_sqv1    = sub_sqv1      , det_thl1    = det_thl1     , &
            det_sqv1        = det_sqv1      ,                                                            &
            maxwidth        = maxwidth1     , maxmf       = maxmf1        , ztop_plume  = ztop_plume1  , &
            excess_h        = excess_h1     , excess_q    = excess_q1     , maxwidth_dd = maxwidth_dd1 , &
            maxmf_dd        = maxmf_dd1     , maxtkeprod  =maxtkeprod1    , cldtop_cooling=cldtop_cooling1,&
            ent_eff         = ent_eff1      ,                                                            &
            flag_qc         = f_qc          , flag_qi     = f_qi          , flag_qs     = f_qs         , &
            flag_ozone      = f_qoz         , flag_qnc    = f_qnc         , flag_qni    = f_qni        , &
            flag_qnwfa      = f_qnwfa       , flag_qnifa  = f_qnifa       , flag_qnbca  = f_qnbca      , &
            pattern_spp_pbl1= pattern_spp1  ,                                                            &
            mix_chem        = mix_chem1     , enh_mix     = enh_mix1      , nchem       = nchem1       , &
            ndvel           = ndvel1        , chem1       = chem1         , emis_ant_no = emisant_no1  , &
            frp             = frp1          , vdep        = vd1           , settle1     = settle1      , &
            nscalars        = nscalars      , scalars     = scalars       ,                              &
            bl_mynn_tkeadvect  = bl_mynn_tkeadvect    , &
            tke_budget         = tke_budget           , &
            bl_mynn_cloudpdf   = bl_mynn_cloudpdf     , &
            bl_mynn_mixlength  = bl_mynn_mixlength    , &
            bl_mynn_closure    = bl_mynn_closure      , &
            bl_mynn_edmf       = bl_mynn_edmf         , &
            bl_mynn_edmf_dd    = bl_mynn_edmf_dd      , &
            bl_mynn_edmf_mom   = bl_mynn_edmf_mom     , &
            bl_mynn_edmf_tke   = bl_mynn_edmf_tke     , &
            bl_mynn_mixscalars = bl_mynn_mixscalars   , &
            bl_mynn_mixaerosols= bl_mynn_mixaerosols  , &
            bl_mynn_mixnumcon  = bl_mynn_mixnumcon    , &
            bl_mynn_output     = bl_mynn_output       , &
            bl_mynn_cloudmix   = bl_mynn_cloudmix     , &
            bl_mynn_mixqt      = bl_mynn_mixqt        , &
            bl_mynn_ess        = bl_mynn_ess          , &
            spp_pbl            = spp_pbl              , &
            kts = kts , kte = kte , errmsg = errmsg , errflg = errflg )

    if (debug) then
       print*,"In mynnedmf driver, after call to mynnedmf"
    endif

    !--- conversion of tendencies in terms of specific contents to in terms of mixing ratios:
    if (dry_mixing_ratio) then
       call  mynnedmf_post_run(kte,f_qc,f_qi,f_qs,delt,qv1,qc1,qi1,qs1,rqvblten1,rqcblten1, &
                              rqiblten1,rqsblten1,errmsg,errflg)
    endif

    !--- inout arguments:
    pblh(i,j)  = pblh1
    kpbl(i,j)  = kpbl1

   if (dry_mixing_ratio) then
      do k=kts,kte
         qc_bl(i,k,j)     = qc_bl1(k)/(one - sqv1(k))
         qi_bl(i,k,j)     = qi_bl1(k)/(one - sqv1(k))
         cldfra_bl(i,k,j) = cldfra_bl1(k)
      enddo
   else
      do k=kts,kte
         qc_bl(i,k,j)     = qc_bl1(k)
         qi_bl(i,k,j)     = qi_bl1(k)
         cldfra_bl(i,k,j) = cldfra_bl1(k)
      enddo
   endif

    do k = kts,kte
       el_pbl(i,k,j)  = el_pbl1(k)
       qke(i,k,j)     = qke1(k)
       qke_adv(i,k,j) = qke_adv1(k)
       cov(i,k,j)     = cov1(k)
       tsq(i,k,j)     = tsq1(k)
       qsq(i,k,j)     = qsq1(k)
       sh3d(i,k,j)    = sh1(k)
       sm3d(i,k,j)    = sm1(k)
       exch_h(i,k,j)  = kh1(k)
       exch_m(i,k,j)  = km1(k)
    enddo

    !--- inout tendencies:
    do k = kts,kte
       rublten(i,k,j)    = rublten1(k) 
       rvblten(i,k,j)    = rvblten1(k) 
       rthblten(i,k,j)   = rthblten1(k) 
    enddo

    if (present(rqvblten)) then
       do k=kts,kte
          rqvblten(i,k,j) = rqvblten1(k)
       enddo
    end if
    if(f_qc .and. present(rqcblten)) then
       do k = kts,kte
          rqcblten(i,k,j) = rqcblten1(k) 
       enddo
    endif
    if(f_qi .and. present(rqiblten)) then
       do k = kts,kte
          rqiblten(i,k,j) = rqiblten1(k) 
       enddo
    endif
    if(f_qs .and. present(rqsblten)) then
       do k = kts,kte
          rqsblten(i,k,j) = rqsblten1(k)
       enddo
    endif
    if(f_qoz .and. present(rqozblten)) then
       do k = kts,kte
          rqozblten(i,k,j) = rqozblten1(k) 
       enddo
    endif
    if(f_qnc .and. present(rqncblten)) then
       do k = kts,kte
          rqncblten(i,k,j) = rqncblten1(k) 
       enddo
    endif
    if(f_qni .and. present(rqniblten)) then
       do k = kts,kte
          rqniblten(i,k,j) = rqniblten1(k) 
       enddo
    endif
    if(f_qnifa .and. present(rqnifablten)) then
       do k = kts,kte
          rqnifablten(i,k,j) = rqnifablten1(k) 
       enddo
    endif
    if(f_qnwfa .and. present(rqnwfablten)) then
       do k = kts,kte
          rqnwfablten(i,k,j) = rqnwfablten1(k) 
       enddo
    endif
    if(f_qnbca .and. present(rqnbcablten)) then
       do k = kts,kte
          rqnbcablten(i,k,j) = rqnbcablten1(k) 
       enddo
    endif

    if (bl_mynn_output > 0) then
      do k = kts,kte
         edmf_a(i,k,j)   = edmf_a1(k)
         edmf_w(i,k,j)   = edmf_w1(k)
         edmf_qt(i,k,j)  = edmf_qt1(k)
         edmf_thl(i,k,j) = edmf_thl1(k)
         edmf_ent(i,k,j) = edmf_ent1(k)
         edmf_qc(i,k,j)  = edmf_qc1(k)
         sub_thl(i,k,j)  = sub_thl1(k)
         sub_sqv(i,k,j)  = sub_sqv1(k)
         det_thl(i,k,j)  = det_thl1(k)
         det_sqv(i,k,j)  = det_sqv1(k)
      enddo
    endif

    !--- output arguments:
    if (bl_mynn_edmf > 0) then
      maxwidth(i,j)      = maxwidth1
      maxmf(i,j)         = maxmf1
      ztop_plume(i,j)    = ztop_plume1
      excess_h(i,j)      = excess_h1
      excess_q(i,j)      = excess_q1
      maxwidth_dd(i,j)   = maxwidth_dd1
      maxmf_dd(i,j)      = maxmf_dd1
      maxtkeprod(i,j)    = maxtkeprod1
      cldtop_cooling(i,j)= cldtop_cooling1
      ent_eff(i,j)       = ent_eff1
    endif

    if( (tke_budget .eq. 1) .and. present(qwt)   .and. present(qbuoy) .and. present(qshear) .and. &
       present(qdiss) .and. present(dqke)) then
       do k = kts,kte
          dqke(i,k,j)   = dqke1(k)
          qwt(i,k,j)    = qwt1(k)
          qshear(i,k,j) = qshear1(k)
          qbuoy(i,k,j)  = qbuoy1(k)
          qdiss(i,k,j)  = qdiss1(k)
       enddo
    endif

    if (mix_chem1 .and. present(chem3d)) then
       do nc = 1,nchem
          do k = kts,kte
             chem3d(i,k,j,nc) = max(1.e-12, chem1(k,nc))
          enddo
       enddo
    endif
    deallocate(vd1)
    deallocate(chem1)
    deallocate(settle1)

    !--- calculating MYNN-EDMF diagnostics:
    if (debug) then
       write(0,*)"bl_mynn_diags=", bl_mynn_diags
       write(0,*)"In mynnedmf driver, just before call to mynnedmf_diags"
    endif

    if (bl_mynn_diags >= 1) then
       call mynnedmf_diags (&
               kts  = kts , kte    = kte    , delp1      = delp1     , dz1  =  dz1  , zw1 = zw1         ,&
               zagl1=zagl1, u1     = u1     ,                                                            &
               v1   = v1  , tk1    = tk1    , qc1        = qc1       , qi1  =  qi1                      ,&
               qs1  = qs1 , qc_bl1 = qc_bl1 , qi_bl1     = qi_bl1    , cldfra_bl1 = cldfra_bl1          ,&        
               rho1 = rho1, xland1 = xland1 , pblh1      = pblh1     ,                                   &
               ! diagnostic outputs
               lwp1     = lwp1   , iwp1    = iwp1   , swp1       = swp1      , cldceil1 = cldceil1      ,&
               wspd101  = wspd101, wspd801 = wspd801, wspd1601   = wspd1601  , maxcldfra1 = maxcldfra1  ,&
               maxcldfra_pbl1 = maxcldfra_pbl1      , bl_mynn_diags = bl_mynn_diags )

       ! collect diagnostic output
       lwp(i,j)     = lwp1
       iwp(i,j)     = iwp1
       swp(i,j)     = swp1
       cldceil(i,j) = cldceil1

       if (bl_mynn_diags >= 2) then
          wspd10(i,j)  = wspd101
          wspd80(i,j)  = wspd801
          wspd160(i,j) = wspd1601

          maxcldfra(i,j)     = maxcldfra1
          maxcldfra_pbl(i,j) = maxcldfra_pbl1
       endif

    endif

 enddo !i
 enddo !j

 if (tke_budget .eq. 1) then
   deallocate(dqke1     )
   deallocate(qwt1      )
   deallocate(qshear1   )
   deallocate(qbuoy1    )
   deallocate(qdiss1    )
 endif

 if (debug) then
   print*,"In mynnedmf_driver, at end"
 endif
 
 end subroutine mynnedmf_driver

! ==================================================================
 subroutine moisture_check2(kte,  delt,    dp, exner, &
                             qv,    qc,    qi,    th  )
  !
  ! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  ! force them to be larger than minimum value by (1) condensating 
  ! water vapor into liquid or ice, and (2) by transporting water vapor 
  ! from the very lower layer.
  ! 
  ! We then update the final state variables and tendencies associated
  ! with this correction. If any condensation happens, update theta/temperature too.
  ! Note that (qv,qc,qi,th) are the final state variables after
  ! applying corresponding input tendencies and corrective tendencies.

    use module_bl_mynnedmf_common, only: kind_phys,xlvcp,xlscp,zero,two,p5

    implicit none
    integer,         intent(in)     :: kte
    real(kind_phys), intent(in)     :: delt
    real(kind_phys), dimension(kte), intent(in)     :: dp
    real(kind_phys), dimension(kte), intent(in)     :: exner
    real(kind_phys), dimension(kte), intent(inout)  :: qv, qc, qi, th
    integer   k
    real(kind_phys) ::  dqc2, dqi2, dqv2, sum, aa, dum
    real(kind_phys), parameter :: qvmin1= 1e-8,    & !min at k=1
                                  qvmin = 1e-20,   & !min above k=1
                                  qcmin = 0.0,     &
                                  qimin = 0.0

    do k = kte, 1, -1  ! From the top to the surface
       dqc2 = max(zero, qcmin-qc(k)) !qc deficit (>=0)
       dqi2 = max(zero, qimin-qi(k)) !qi deficit (>=0)

       !update species
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qv(k)  = qv(k)  -  dqc2 - dqi2
       !for theta
       th(k)  = th(k)  +  xlvcp/exner(k)*dqc2 + &
                          xlscp/exner(k)*dqi2
       !for temperature
       !th(k)  = th(k)  +  xlvcp*dqc2 + &
       !                   xlscp*dqi2

       !then fix qv if lending qv made it negative
       if (k .eq. 1) then
          dqv2   = max(zero, qvmin1-qv(k)) !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k)  = max(qv(k),qvmin1)
          dqv2   = zero
       else
          dqv2   = max(zero, qvmin-qv(k))  !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k-1)= qv(k-1)  - dqv2*dp(k)/dp(k-1)
          qv(k)  = max(qv(k),qvmin)
       endif
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
    end do

    ! Extra moisture used to satisfy 'qv(1)>=qvmin' is proportionally
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture.
    if( dqv2 .gt. 1.e-20 ) then
        sum = zero
        do k = 1, kte
           if( qv(k) .gt. two*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20_kind_phys,sum)
        if( aa .lt. p5 ) then
            do k = 1, kte
               if( qv(k) .gt. two*qvmin ) then
                   dum    = aa*qv(k)
                   qv(k)  = qv(k) - dum
               endif
            enddo
        else
        ! For testing purposes only (not yet found in any output):
        !    write(*,*) 'Full moisture conservation is impossible'
        endif
    endif

    return

 end subroutine moisture_check2
!=================================================================================================================
!>\section arg_table_mynnedmf_pre_init
!!\html\include mynnedmf_pre_init.html
!!
 subroutine mynnedmf_pre_init(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_init

!=================================================================================================================
!>\section arg_table_mynnedmf_pre_finalize
!!\html\include mynnedmf_pre_finalize.html
!!
 subroutine mynnedmf_pre_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_finalize

!=================================================================================================================
!>\section arg_table_mynnedmf_pre_run
!!\html\include mynnedmf_pre_run.html
!!
 subroutine mynnedmf_pre_run(kte,f_qc,f_qi,f_qs,qv,qc,qi,qs,sqv,sqc,sqi,sqs,errmsg,errflg)
!=================================================================================================================
 use module_bl_mynnedmf_common,only: kind_phys,zero,one
   
!--- input arguments:
 logical,intent(in):: &
    f_qc,      &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi,      &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs        ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: kte

 real(kind_phys),intent(in),dimension(1:kte):: &
    qv,        &!
    qc,        &!
    qi,        &!
    qs          !


!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

 real(kind_phys),intent(out),dimension(1:kte):: &
    sqv,       &!
    sqc,       &!
    sqi,       &!
    sqs         !


!--- local variables:
 integer:: k
 integer,parameter::kts=1
!-----------------------------------------------------------------------------------------------------------------

!--- initialization:
 do k = kts,kte
    sqc(k) = zero
    sqi(k) = zero
 enddo

!--- conversion from water vapor mixing ratio to specific humidity:
 do k = kts,kte
    sqv(k) = qv(k)/(one+qv(k))
 enddo

!--- conversion from cloud liquid water,cloud ice,and snow mixing ratios to specific contents:
 if(f_qc) then
    do k = kts,kte
       sqc(k) = qc(k)/(one+qv(k))
    enddo
 endif
 if(f_qi) then
    do k = kts,kte
       sqi(k) = qi(k)/(one+qv(k))
    enddo
 endif
 if(f_qs) then
    do k = kts,kte
       sqs(k) = qs(k)/(one+qv(k))
    enddo
 endif

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_run
!=================================================================================================================

!=================================================================================================================
!>\section arg_table_mynnedmf_post_init
!!\html\include mynnedmf_post_init.html
!!
 subroutine mynnedmf_post_init(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_post_init

!=================================================================================================================
!>\section arg_table_mynnedmf_post_finalize
!!\html\include mynnedmf_post_finalize.html
!!
 subroutine mynnedmf_post_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_post_finalize

!=================================================================================================================
!>\section arg_table_mynnedmf_post_run
!!\html\include mynnedmf_post_run.html
!!
 subroutine mynnedmf_post_run(kte,f_qc,f_qi,f_qs,delt,qv,qc,qi,qs,dqv,dqc,dqi,dqs,errmsg,errflg)
!=================================================================================================================
 use module_bl_mynnedmf_common,only: kind_phys,one
!--- input arguments:
 logical,intent(in):: &
    f_qc, &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi, &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs   ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: kte

 real(kind_phys),intent(in):: &
    delt   !

 real(kind_phys),intent(in),dimension(1:kte):: &
    qv,   &!
    qc,   &!
    qi,   &!
    qs     !


!--- inout arguments:
 real(kind_phys),intent(inout),dimension(1:kte):: &
    dqv,  &!
    dqc,  &!
    dqi,  &!
    dqs    !


!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg


!--- local variables:
 integer:: k
 integer,parameter::kts=1
 real(kind_phys):: rq,sq,tem
 real(kind_phys),dimension(1:kte):: sqv,sqc,sqi,sqs

!-----------------------------------------------------------------------------------------------------------------

!--- initialization:
 do k = kts,kte
    sq = qv(k)/(one+qv(k))     !conversion of qv at time-step n from mixing ratio to specific humidity.
    sqv(k) = sq + dqv(k)*delt  !calculation of specific humidity at time-step n+1.
    rq = sqv(k)/(one-sqv(k))   !conversion of qv at time-step n+1 from specific humidity to mixing ratio.
    dqv(k) = (rq - qv(k))/delt !calculation of the tendency.
 enddo

 if(f_qc) then
    do k = kts,kte
       sq = qc(k)/(one+qv(k))
       sqc(k) = sq + dqc(k)*delt
       rq  = sqc(k)/(one-sqv(k))
       dqc(k) = (rq - qc(k))/delt
    enddo
 endif

 if(f_qi) then
    do k = kts,kte
       sq = qi(k)/(one+qv(k))
       sqi(k) = sq + dqi(k)*delt
       rq = sqi(k)/(one-sqv(k))
       dqi(k) = (rq - qi(k))/delt
    enddo
 endif

 if(f_qs) then
    do k = kts,kte
       sq = qs(k)/(one+qv(k))
       sqs(k) = sq + dqs(k)*delt
       rq = sqs(k)/(one-sqv(k))
       dqs(k) = (rq - qs(k))/delt
    enddo
 endif

!--- output error flag and message:
 errmsg = " "
 errflg = 0

 end subroutine mynnedmf_post_run

!=================================================================================================================
 end module module_bl_mynnedmf_driver
!=================================================================================================================

