!> \file module_bl_mynnedmf_driver.F90
!!  This serves as the unified interface between the WRF/MPAS PBL drivers and the MYNN-EDMF
!!  turbulence scheme (in module_bl_mynnedmf.F90).

!>\ingroup gsd_mynn_edmf
!> The following references best describe the code within
!!    Olson et al. (2026, NOAA Technical Memorandum)
!!    Nakanishi and Niino (2009) \cite NAKANISHI_2009
!=================================================================================================================
 module module_bl_mynnedmf_driver

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

   implicit none

   logical,intent(in) :: allowed_to_read,restart

   integer,intent(in) :: ids,ide,jds,jde,kds,kde,  &
        &                ims,ime,jms,jme,kms,kme,  &
        &                its,ite,jts,jte,kts,kte

   real(kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout) :: &
        &rublten,rvblten,rthblten,rqvblten,        &
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
! \brief This subroutine (1) performs preparational work (conversions of units and assemble column arrays),
!!  (2) runs the mynnedmf scheme, and then (3) performs post-mynnedmf work (conversions and updates 3d arrays).
!> \section arg_table_mynnedmf_driver Argument Table
!! \htmlinclude mynnedmf_driver.html
!!
 SUBROUTINE mynnedmf_driver(           &
                  ids               , ide               , jds                , jde                , &
                  kds               , kde               , ims                , ime                , &
                  jms               , jme               , kms                , kme                , &
                  its               , ite               , jts                , jte                , &
                  kts               , kte               , flag_qc            , flag_qi            , &
                  flag_qs           , flag_qnc          , flag_qni           , flag_oz            , &
                  flag_qnifa        , flag_qnwfa        , flag_qnbca         , initflag           , &
                  restart           , cycling           , delt               ,                      &
                  dx                , xland             , ps                 , ts                 , &
                  qsfc              , ust               , ch                 , hfx                , &
                  qfx               , wspd              , znt                ,                      &
                  uoce              , voce              ,                                           &
                  !3d input
                  dz                , u                 , v                  , w                  , &
                  th                , t3d               , p                  , exner              , &
                  rho               , qv                , qc                 , qi                 , &
                  qs                , qnc               , qni                , qnifa              , &
                  qnwfa             , qnbca             , qoz                , rthraten           , &
                  !3d output
                  cldfra_bl         , qc_bl             , qi_bl              ,                      &
                  qke               , qke_adv           , el_pbl             , sh3d               , &
                  tsq               , qsq               , cov                , sm3d               , &
                  !output tendencies
                  rublten           , rvblten           , rthblten           ,                      &
                  rqvblten          , rqcblten          , rqiblten           , rqsblten           , &
                  rqncblten         , rqniblten         , rqnifablten        , rqnwfablten        , &
                  rqnbcablten       , rqozblten         ,                                           &
                  !2d output
                  pblh              , kpbl              , maxwidth           ,                      &
                  maxmf             , ztop_plume        , excess_h           , excess_q           , &
                  !optional 3d output
                  edmf_a            , edmf_w            ,                                           &
                  edmf_qt           , edmf_thl          , edmf_ent           , edmf_qc            , &
                  sub_thl           , sub_sqv           , det_thl            , det_sqv            , &
                  exch_h            , exch_m            , dqke               , qwt                , &
                  qshear            , qbuoy             , qdiss              ,                      &
                  !configuration options (+spp array)
                  spp_pbl           , pattern_spp       , icloud_bl          ,                      &
                  bl_mynn_tkeadvect , tke_budget        , bl_mynn_cloudpdf   , bl_mynn_mixlength  , &
                  bl_mynn_closure   , bl_mynn_edmf      , bl_mynn_edmf_mom   , bl_mynn_edmf_tke   , &
                  bl_mynn_output    , bl_mynn_mixscalars, bl_mynn_mixaerosols, bl_mynn_mixnumcon  , &
                  bl_mynn_cloudmix  , bl_mynn_mixqt     , bl_mynn_edmf_dd    , bl_mynn_ess        , &
                  !smoke/dust
                  mix_chem          , nchem             , ndvel              , enh_mix            , &
                  chem3d            , settle3d          , vd3d               ,                      &
                  frp_mean          , emis_ant_no       ,                                           &
                  !ccpp error handling
                  errmsg            , errflg                                                        &
               )

!=================================================================================================================
!--- input arguments:
 logical,intent(in):: &
    flag_qc,               &     ! if true,the physics package includes the cloud liquid water mixing ratio.
    flag_qi,               &     ! if true,the physics package includes the cloud ice mixing ratio.
    flag_qs,               &     ! if true,the physics package includes the snow mixing ratio.
    flag_qnc,              &     ! if true,the physics package includes the cloud liquid water number concentration.
    flag_qni,              &     ! if true,the physics package includes the cloud ice number concentration.
    flag_qnifa,            &     ! if true,the physics package includes the "ice-friendly" aerosol number concentration.
    flag_qnwfa,            &     ! if true,the physics package includes the "water-friendly" aerosol number concentration.
    flag_qnbca,            &     ! if true,the physics package includes the number concentration of black carbon.
    flag_oz                      ! if true,the physics package includes the number concentration of ozone.
    
 logical, intent(in) ::    &
    bl_mynn_tkeadvect,     &
    cycling,               &
    restart
 
 integer,intent(in):: &
    ids,ide,jds,jde,kds,kde, &
    ims,ime,jms,jme,kms,kme, &
    its,ite,jts,jte,kts,kte

 integer, intent(in) ::    &
    bl_mynn_cloudpdf,      &
    bl_mynn_mixlength,     &
    bl_mynn_edmf,          &
    bl_mynn_edmf_dd,       &
    bl_mynn_edmf_mom,      &
    bl_mynn_edmf_tke,      &
    bl_mynn_cloudmix,      &
    bl_mynn_mixqt,         &
    bl_mynn_output,        &
    bl_mynn_mixscalars,    &
    bl_mynn_mixaerosols,   &
    bl_mynn_mixnumcon,     &
    bl_mynn_ess,           &
    tke_budget

 integer,intent(in):: &
    initflag,              &!
    icloud_bl,             &!
    spp_pbl
 
 real(kind_phys),intent(in):: &
    bl_mynn_closure

 real(kind_phys), intent(in) :: delt

!required 2D input
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

!required 3D input
 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz,          &!
    u,           &!
    w,           &!
    v,           &!
    th,          &!
    t3d,         &!
    p,           &!
    exner,       &!
    rho,         &!
    rthraten      !

!optional 3D input
 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    qv,          &!
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

 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
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

!--- output arguments imposed by ccpp:
 character(len=*),intent(out) :: &
    errmsg        ! output error message (-).

 integer,intent(out):: &
    errflg        ! output error flag (-).

 real(kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    maxwidth,    &!
    maxmf,       &!
    ztop_plume,  &!
    excess_h,    &!
    excess_q

 real(kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    exch_h,      &!
    exch_m        !

 real(kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    dqke,        &!
    qwt,         &!
    qshear,      &!
    qbuoy,       &!
    qdiss         !
 
 real(kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: &
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

!--smoke/chem arrays
 logical,intent(in):: mix_chem
 integer,intent(in):: nchem,ndvel

 real(kind_phys),intent(in),   dimension(ims:ime,jms:jme),      optional:: frp_mean,emis_ant_no
 real(kind_phys),intent(in),   dimension(ims:ime,jms:jme,ndvel),optional:: vd3d
 real(kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme,nchem),optional:: chem3d,settle3d
 logical,intent(inout),optional:: enh_mix

 real(kind_phys),dimension(ndvel):: vd1
 real(kind_phys),dimension(kts:kte,nchem):: chem1,settle1

!--Generic scalar array support (not yet connected to all model systems, but any new scalars that need to be mixed
!  (locally and nonlocally) can be added to this generic array and will be mixed when bl_mynn_mixscalars=1.
 integer, parameter :: nscalars=1
 real(kind_phys),dimension(kts:kte,nscalars):: scalars

!--Local (1d arrays and single-point variables)
 integer :: i, j, k, itf, jtf, n, nc, nd

 integer:: kpbl1

 real(kind_phys):: &
    dx1,xland1,ps1,ts1,qsfc1,ust1,ch1,hfx1,qfx1, &
    wspd1,uoce1,voce1,znt1

 real(kind_phys),dimension(kts:kte):: &
    dz1,u1,v1,th1,tk1,p1,ex1,rho1,qv1,rthraten1

 real(kind_phys),dimension(kts:kme):: &
    w1

 real(kind_phys),dimension(kts:kte):: &
    qc1,qi1,qs1,qnc1,qni1,qnifa1,qnwfa1,qnbca1,qoz1

 real(kind_phys),dimension(kts:kte):: &
    pattern_spp1

 real(kind_phys):: &
    pblh1

 real(kind_phys),dimension(kts:kte):: &
    cldfra_bl1,qc_bl1,qi_bl1,el1,qke1,qke_adv1,cov1,qsq1,tsq1,sh1,sm1,kh1,km1

 real(kind_phys),dimension(kts:kte):: &
    du1,dv1,dth1,dqv1,dqc1,dqi1,dqs1, &
    dqnc1,dqni1,dqnifa1,dqnwfa1,dqnbca1,doz1

 real(kind_phys),dimension(kts:kte):: &
    edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1, &
    sub_thl1,sub_sqv1,det_thl1,det_sqv1

 real(kind_phys):: &
    maxwidth1,maxmf1,ztop_plume1,excess_h1,excess_q1

 real(kind_phys),dimension(:),allocatable:: &
    dqke1,qwt1,qshear1,qbuoy1,qdiss1

 real(kind_phys),dimension(kts:kte):: &
    sqv1,sqc1,sqi1,sqs1

 real(kind_phys):: frp1,emis1

 real(kind_phys), dimension(kts:kte):: kzero

 !for debug printing, set to true
 logical, parameter:: debug = .false.

!-----------------------------------------------------------------------------------------------------------------
 if (debug) then
    write(0,*)"=============================================="
    write(0,*)"in mynn-edmf driver..."
    write(0,*)"initflag=",initflag," restart =",restart
 endif

 errmsg = " "
 errflg = 0
   
 jtf=MIN0(JTE,JDE-1)
 itf=MIN0(ITE,IDE-1)

 !For now, initialized bogus array
 !ozone            =zero
 !rO3blten         =zero
 kzero            =zero
 !initialize subgrid clouds:
 qc_bl1           =zero
 qi_bl1           =zero
 cldfra_bl1       =zero
 !spp
 pattern_spp1     =zero
 !turbulence properties
 qke1             =zero
 qke_adv1         =zero
 el1              =zero
 sh1              =zero
 sm1              =zero
 kh1              =zero
 km1              =zero
 tsq1             =zero
 qsq1             =zero
 cov1             =zero
 !tke budget (optional arrays)
 if (tke_budget .eq. 1) then
    allocate(dqke1(kts:kte),    source=zero)
    allocate(qwt1(kts:kte),     source=zero)
    allocate(qshear1(kts:kte),  source=zero)
    allocate(qbuoy1(kts:kte),   source=zero)
    allocate(qdiss1(kts:kte),   source=zero)
 endif
 !1d mass-flux arrays - most are used in the scheme, so no toptional
 edmf_a1          =zero
 edmf_w1          =zero
 edmf_qt1         =zero
 edmf_thl1        =zero
 edmf_ent1        =zero
 edmf_qc1         =zero
 sub_thl1         =zero
 sub_sqv1         =zero
 det_thl1         =zero
 det_sqv1         =zero
 !moist species
 qv1              =zero
 qc1              =zero
 qi1              =zero
 qs1              =zero
 qnc1             =zero
 qni1             =zero
 qnwfa1           =zero
 qnifa1           =zero
 qnbca1           =zero
 qoz1             =zero
 !1d (non-optional) tendencies
 du1              =zero
 dv1              =zero
 dth1             =zero
 dqv1             =zero
 dqc1             =zero
 dqi1             =zero
 dqs1             =zero
 dqni1            =zero
 dqnc1            =zero
 dqnwfa1          =zero
 dqnifa1          =zero
 dqnbca1          =zero
 doz1             =zero

 !---------------------------------------
 !Begin looping in the i- and j-direction
 !---------------------------------------
 do j = jts, jte !jtf
   do i = its, ite !itf
      !3d variables
      do k=kts,kte
         u1(k)       = u(i,k,j)
         v1(k)       = v(i,k,j)
         w1(k)       = w(i,k,j) !not iitialized at kte+1
         th1(k)      = th(i,k,j)
         p1(k)       = p(i,k,j)
         ex1(k)      = exner(i,k,j)
         rho1(k)     = rho(i,k,j)
         tk1(k)      = t3d(i,k,j)
         dz1(k)      = dz(i,k,j)
         rthraten1(k)= rthraten(i,k,j)
      enddo
      !2d variables
      dx1            = dx(i,j)
      xland1         = xland(i,j)
      ts1            = ts(i,j)
      qsfc1          = qsfc(i,j)
      ps1            = ps(i,j)
      ust1           = ust(i,j)
      ch1            = ch(i,j)
      wspd1          = wspd(i,j)
      uoce1          = uoce(i,j)
      voce1          = voce(i,j)
      znt1           = znt(i,j)
      !output
      pblh1          = pblh(i,j)
      kpbl1          = kpbl(i,j)
      if (bl_mynn_edmf > 0) then
         maxwidth1      = maxwidth(i,j)
         maxmf1         = maxmf(i,j)
         ztop_plume1    = ztop_plume(i,j)
         excess_h1      = excess_h(i,j)
         excess_q1      = excess_q(i,j)
      endif
      !check for unearthly incoming surface fluxes. These threshold are only surpassed
      !when something unphysical is happening. If these limits are being surpassed,
      !conservation is already questionable and the simulation is heading off the rails. 
      !Try to curb the consequences of this behavior by imposing liberal limits on
      !the incoming fluxes:
      hfx1 = hfx(i,j)
      if (hfx1 > 1200.) then
         !print*,"hfx at i=",i," j=",j,"is unrealistic:",hfx1
         hfx1 = 1200._kind_phys
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
      
      !spp input
      if (spp_pbl > 0) then
         do k=kts,kte
            pattern_spp1(k) = pattern_spp(i,k,j)
         enddo
      endif

      !when NOT cold-starting on the first time step, update input
      if (initflag .eq. 0 .or. restart) THEN
         !update sgs cloud info.
         if (icloud_bl > 0) then
            do k=kts,kte
               qc_bl1(k)     = qc_bl(i,k,j)
               qi_bl1(k)     = qi_bl(i,k,j)
               cldfra_bl1(k) = cldfra_bl(i,k,j)
            enddo
         endif

         !turbulennce variables
         do k=kts,kte
            qke1(k) = qke(i,k,j)
            qsq1(k) = qsq(i,k,j)
            tsq1(k) = tsq(i,k,j)
            cov1(k) = cov(i,k,j)
            sh1(k)  = sh3d(i,k,j)
            sm1(k)  = sm3d(i,k,j)
            kh1(k)  = exch_h(i,k,j)
            km1(k)  = exch_m(i,k,j)
            el1(k)  = el_pbl(i,k,j)
         enddo
         if (bl_mynn_tkeadvect) then
            qke_adv1(kts:kte) = qke_adv(i,kts:kte,j)
         else
            qke_adv1(kts:kte) = qke(i,kts:kte,j)
         endif
      endif
      
      !intialize moist species; guard against incoming negatives
      do k=kts,kte
         qv1(k) = max(1e-10_kind_phys, qv(i,k,j))
      enddo
      if (flag_qc .and. present(qc)) then
         do k=kts,kte
            qc1(k) = max(zero, qc(i,k,j))
         enddo
      endif
      if (flag_qi .and. present(qi)) then
         do k=kts,kte
            qi1(k) = max(zero, qi(i,k,j))
         enddo
      endif
      if (flag_qs .and. present(qs)) then
         do k=kts,kte
            qs1(k) = max(zero, qs(i,k,j))
         enddo
      endif
      if (flag_qnc .and. present(qnc)) then
         do k=kts,kte
            qnc1(k) = qnc(i,k,j)
         enddo
      endif
      if (flag_qni .and. present(qni)) then
         do k=kts,kte
            qni1(k) = qni(i,k,j)
         enddo
      endif
      if (flag_qnwfa .and. present(qnwfa)) then
         do k=kts,kte
            qnwfa1(k) = qnwfa(i,k,j)
         enddo
      endif
      if (flag_qnifa .and. present(qnifa)) then
         do k=kts,kte
            qnifa1(k) = qnifa(i,k,j)
         enddo
      endif
      if (flag_qnbca .and. present(qnbca)) then
         do k=kts,kte
            qnbca1(k) = qnbca(i,k,j)
         enddo
      endif
      if (flag_oz .and. present(qoz)) then
         do k=kts,kte
            qoz1(k) = qoz(i,k,j)
         enddo
      endif

      !Smoke/dust
      if (present(chem3d).and. present(settle3d) .and. present(vd3d) .and. &
         present(frp_mean) .and. present(emis_ant_no)) then
         do nc = 1,nchem
            do k = kts,kte
               chem1(k,nc)   = chem3d(i,k,j,nc)
               settle1(k,nc) = settle3d(i,k,j,nc)
            enddo
         enddo
         do nd = 1,ndvel
            vd1(nd)  = vd3d(i,j,nd)
         enddo
         frp1        = frp_mean(i,j)
         emis1       = emis_ant_no(i,j)
      else
         chem1       = zero
         settle1     = zero
         vd1         = zero
         frp1        = zero
         emis1       = zero
      endif
      
      !generic scalar array support
      scalars     = zero

      !In WRF/MPAS, mixing ratio is incoming; convert to specific contents:
      call mynnedmf_pre_run(kte    , flag_qc , flag_qi , flag_qs ,   &
                            qv1    , qc1     , qi1     , qs1     ,   &
                            sqv1   , sqc1    , sqi1    , sqs1    ,   &
                            errmsg , errflg                          )

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
            pres1           = p1            , ex1         = ex1           , rho1        = rho1         , &
            tk1             = tk1           , xland       = xland1        , ts          = ts1          , &
            qsfc            = qsfc1         , ps          = ps1           , ust         = ust1         , &
            ch              = ch1           , hfx         = hfx1          , qfx         = qfx1         , &
            wspd            = wspd1         , uoce        = uoce1         , voce        = voce1        , &
            qke1            = qke1          , qke_adv1    = qke_adv1      ,                              &
            tsq1            = tsq1          , qsq1        = qsq1          , cov1        = cov1         , &
            rthraten1       = rthraten1     , du1         = du1           , dv1         = dv1          , &
            dth1            = dth1          , dqv1        = dqv1          , dqc1        = dqc1         , &
            dqi1            = dqi1          , dqs1        = kzero         , dqnc1       = dqnc1        , &
            dqni1           = dqni1         , dqnwfa1     = dqnwfa1       , dqnifa1     = dqnifa1      , &
            dqnbca1         = dqnbca1       , dozone1     = doz1          , kh1         = kh1          , &
            km1             = km1           , pblh        = pblh1         , kpbl        = kpbl1        , &
            el1             = el1           , dqke1       = dqke1         , qwt1        = qwt1         , &
            qshear1         = qshear1       , qbuoy1      = qbuoy1        , qdiss1      = qdiss1       , &
            sh1             = sh1           , sm1         = sm1           , qc_bl1      = qc_bl1       , &
            qi_bl1          = qi_bl1        , cldfra_bl1  = cldfra_bl1    ,                              &
            edmf_a1         = edmf_a1       , edmf_w1     = edmf_w1       , edmf_qt1    = edmf_qt1     , &
            edmf_thl1       = edmf_thl1     , edmf_ent1   = edmf_ent1     , edmf_qc1    = edmf_qc1     , &
            sub_thl1        = sub_thl1      , sub_sqv1    = sub_sqv1      , det_thl1    = det_thl1     , &
            det_sqv1        = det_sqv1      ,                                                            &
            maxwidth        = maxwidth1     , maxmf       = maxmf1        , ztop_plume  = ztop_plume1  , &
            excess_h        = excess_h1     , excess_q    = excess_q1     ,                              &
            flag_qc         = flag_qc       , flag_qi     = flag_qi       , flag_qs     = flag_qs      , &
            flag_ozone      = flag_oz       , flag_qnc    = flag_qnc      , flag_qni    = flag_qni     , &
            flag_qnwfa      = flag_qnwfa    , flag_qnifa  = flag_qnifa    , flag_qnbca  = flag_qnbca   , &
            pattern_spp_pbl1= pattern_spp1  , scalars     = scalars       , nscalars    = nscalars     , &
            mix_chem        = mix_chem      , enh_mix     = enh_mix       , nchem       = nchem        , &
            ndvel           = ndvel         , chem1       = chem1         , emis_ant_no = emis1        , &
            frp             = frp1          , vdep        = vd1           , settle1     = settle1      , &
            bl_mynn_tkeadvect  = bl_mynn_tkeadvect    , &
            tke_budget         = tke_budget           , &
            bl_mynn_cloudpdf   = bl_mynn_cloudpdf     , &
            bl_mynn_mixlength  = bl_mynn_mixlength    , &
            closure            = bl_mynn_closure      , &
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
            icloud_bl          = icloud_bl            , &
            spp_pbl            = spp_pbl              , &
            kts = kts , kte = kte , errmsg = errmsg , errflg = errflg )

      if (debug) then
         print*,"In mynnedmf driver, after call to mynnedmf"
      endif

      !--- conversion of tendencies in terms of specific contents to mixing ratios:
      call  mynnedmf_post_run(                                        &
                kte      , flag_qc  , flag_qi  , flag_qs , delt     , &
                qv1      , qc1      , qi1      , qs1     , dqv1     , &
                dqc1     , dqi1     , dqs1     , errmsg  , errflg     )

      ! update turbulence properties output
      do k=kts,kte
         qke(i,k,j)     = qke1(k)
         el_pbl(i,k,j)  = el1(k)
         sh3d(i,k,j)    = sh1(k)
         sm3d(i,k,j)    = sm1(k)
         exch_h(i,k,j)  = kh1(k)
         exch_m(i,k,j)  = km1(k)
         tsq(i,k,j)     = tsq1(k)
         qsq(i,k,j)     = qsq1(k)
         cov(i,k,j)     = cov1(k)
         qke_adv(i,k,j) = qke_adv1(k)
      enddo

      !2d output
      kpbl(i,j)        = kpbl1
      pblh(i,j)        = pblh1
      if (bl_mynn_edmf > 0) then
         maxwidth(i,j)    = maxwidth1
         maxmf(i,j)       = maxmf1
         ztop_plume(i,j)  = ztop_plume1
         excess_h(i,j)    = excess_h1
         excess_q(i,j)	  = excess_q1
      endif

      !- Update 3d tendencies (conversions done above):
      do k=kts,kte
         rublten(i,k,j)  = du1(k)
         rvblten(i,k,j)  = dv1(k)
         rthblten(i,k,j) = dth1(k)
      enddo
      if (present(rqvblten)) then
         do k=kts,kte
            rqvblten(i,k,j) = dqv1(k)
         enddo
      endif
      if (present(rqcblten)) then
         do k=kts,kte
            rqcblten(i,k,j) = dqc1(k)
         enddo
      endif
      if (present(rqiblten)) then
         do k=kts,kte
            rqiblten(i,k,j) = dqi1(k)
         enddo
      endif
      if (present(rqsblten)) then !.false.) then !as of now, there is no rqsblten in wrf
        do k=kts,kte
           rqsblten(i,k,j) = dqs1(k)
        enddo
      endif
      if (present(rqncblten)) then
         do k=kts,kte
            rqncblten(i,k,j) = dqnc1(k)
         enddo
      endif
      if (present(rqniblten)) then
         do k=kts,kte
            rqniblten(i,k,j) = dqni1(k)
         enddo
      endif
      if (present(rqnwfablten)) then
         do k=kts,kte
            rqnwfablten(i,k,j) = dqnwfa1(k)
         enddo
      endif
      if (present(rqnifablten)) then
         do k=kts,kte
            rqnifablten(i,k,j) = dqnifa1(k)
         enddo
      endif
      if (present(rqnbcablten)) then
         do k=kts,kte
            rqnbcablten(i,k,j) = dqnbca1(k)
         enddo
      endif

     !- Collect 3D ouput:
      if (icloud_bl > 0) then
         do k=kts,kte
            qc_bl(i,k,j)     = qc_bl1(k)/(one - sqv1(k))
            qi_bl(i,k,j)     = qi_bl1(k)/(one - sqv1(k))
            cldfra_bl(i,k,j) = cldfra_bl1(k)
         enddo
      endif

      if (tke_budget .eq. 1) then
         do k=kts,kte
            dqke(i,k,j)      = dqke1(k)
            qwt(i,k,j)       = qwt1(k)
            qshear(i,k,j)    = qshear1(k)
            qbuoy(i,k,j)     = qbuoy1(k)
            qdiss(i,k,j)     = qdiss1(k)
         enddo
      endif

      if (bl_mynn_output > 0) then
         do k=kts,kte
            edmf_a(i,k,j)    = edmf_a1(k)
            edmf_w(i,k,j)    = edmf_w1(k)
            edmf_qt(i,k,j)   = edmf_qt1(k)
            edmf_thl(i,k,j)  = edmf_thl1(k)
            edmf_ent(i,k,j)  = edmf_ent1(k)
            edmf_qc(i,k,j)   = edmf_qc1(k)
            sub_thl(i,k,j)   = sub_thl1(k)
            sub_sqv(i,k,j)   = sub_sqv1(k)
            det_thl(i,k,j)   = det_thl1(k)
            det_sqv(i,k,j)   = det_sqv1(k)
         enddo
      endif

      if (mix_chem .and. present(chem3d)) then
         do n = 1,nchem
            do k = kts,kte
               chem3d(i,k,j,n) = max(1.e-12_kind_phys, chem1(k,n))
            enddo
         enddo
      endif

   enddo  !end j-loop
   enddo  !end i-loop

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


!=================================================================================================================
!>\section arg_table_mynnedmf_pre_run
!!\html\include mynnedmf_pre_run.html
!!
 subroutine mynnedmf_pre_run(kte,f_qc,f_qi,f_qs,qv,qc,qi,qs,sqv,sqc,sqi,sqs,errmsg,errflg)
!=================================================================================================================

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

!---  output arguments:
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
!>\section arg_table_mynnedmf_post_run
!!\html\include mynnedmf_post_run.html
!!
 subroutine mynnedmf_post_run(kte,f_qc,f_qi,f_qs,delt,qv,qc,qi,qs,dqv,dqc,dqi,dqs,errmsg,errflg)
!=================================================================================================================

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
 real(kind_phys):: rq,sq
 real(kind_phys),dimension(1:kte):: sqv,sqc,sqi,sqs
!-----------------------------------------------------------------------------------------------------------------
!---  initialization:
 do k = kts,kte
    sq = qv(k)/(one+qv(k))     !conversion of qv at time-step n from mixing ratio to specific humidity.
    sqv(k) = sq + dqv(k)*delt  !calculation of specific humidity at time-step n+1.
    rq = sqv(k)/(one-sqv(k))   !conversion of qv at time-step n+1 from specific humidity to mixing ratio.
    dqv(k) = (rq - qv(k))/delt !calculation of the tendency.
 enddo

 if (f_qc) then
    do k = kts,kte
       sq = qc(k)/(one+qv(k))
       sqc(k) = sq + dqc(k)*delt
       rq  = sqc(k)/(one-sqv(k))
       !rq = sqc(k)*(one+sqv(k))
       dqc(k) = (rq - qc(k))/delt
    enddo
 endif

 if (f_qi) then
    do k = kts,kte
       sq = qi(k)/(one+qv(k))
       sqi(k) = sq + dqi(k)*delt
       rq = sqi(k)/(one-sqv(k))
       !rq = sqi(k)*(one+sqv(k))
       dqi(k) = (rq - qi(k))/delt
    enddo
 endif

 if (f_qs) then
    do k = kts,kte
       sq = qs(k)/(one+qv(k))
       sqs(k) = sq + dqs(k)*delt
       rq = sqs(k)/(one-sqv(k))
       !rq = sqs(k)*(one+sqv(k))
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

