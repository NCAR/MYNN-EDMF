!>\file module_bl_mynnedmf_common.F90
!! Define Model-specific constants/parameters.
!! This module will be used at the initialization stage
!! where all model-specific constants are read and saved into
!! memory. This module is then used again in the MYNN-EDMF. All
!! MYNN-specific constants are declared globally in the main 
!! module (module_bl_mynn) further below:

!>\ingroup gp_mynnedmf
!! This module defines model-specific constants/parameters.
 module module_bl_mynnedmf_common

!------------------------------------------
!
!------------------------------------------

! The following 5-6 lines are the only lines in this file that are not 
! universal for all dycores... Any ideas how to universalize it?
! For MPAS:
! use mpas_kind_types,only: kind_phys => RKIND
! For CCPP:
  use machine,  only : kind_phys

 implicit none
 save

! To be specified from dycore
 real(kind_phys):: cp           !<= 7.*r_d/2. (J/kg/K)
 real(kind_phys):: cpv          !<= 4.*r_v    (J/kg/K) Spec heat H2O gas
 real(kind_phys):: cice         !<= 2106.     (J/kg/K) Spec heat H2O ice
 real(kind_phys):: cliq         !<= 4190.     (J/kg/K) Spec heat H2O liq
 real(kind_phys):: p608         !<= R_v/R_d-1.
 real(kind_phys):: ep_2         !<= R_d/R_v
 real(kind_phys):: grav         !<= accel due to gravity
 real(kind_phys):: karman       !<= von Karman constant
 real(kind_phys):: t0c          !<= temperature of water at freezing, 273.15 K
 real(kind_phys):: rcp          !<= r_d/cp
 real(kind_phys):: r_d          !<= 287.  (J/kg/K) gas const dry air
 real(kind_phys):: r_v          !<= 461.6 (J/kg/K) gas const water
 real(kind_phys):: xlf          !<= 0.35E6 (J/kg) fusion at 0 C
 real(kind_phys):: xlv          !<= 2.50E6 (J/kg) vaporization at 0 C
 real(kind_phys):: xls          !<= 2.85E6 (J/kg) sublimation
 real(kind_phys):: rvovrd       !<= r_v/r_d != 1.608

! Specified locally
 real(kind_phys),parameter:: tref  = 300.0     !< reference temperature (K)
 real(kind_phys),parameter:: TKmin = 253.0     !< for total water conversion, Tripoli and Cotton (1981)
 real(kind_phys),parameter:: p1000mb=100000.0
 real(kind_phys),parameter:: svp1  = 0.6112    !<(kPa)
 real(kind_phys),parameter:: svp2  = 17.67     !<(dimensionless)
 real(kind_phys),parameter:: svp3  = 29.65     !<(K)
 real(kind_phys),parameter:: tice  = 240.0     !<-33 (C), temp at saturation w.r.t. ice
 real(kind_phys),parameter:: tliq  = 269.0     !all hydrometeors are liquid when T > tliq

! To be derived in the init routine
 real(kind_phys):: ep_3         !<= 1.-ep_2 != 0.378
 real(kind_phys):: gtr          !<= grav/tref
 real(kind_phys):: rk           !<= cp/r_d
 real(kind_phys):: tv0          !<= p608*tref
 real(kind_phys):: tv1          !<= (1.+p608)*tref
 real(kind_phys):: xlscp        !<= (xlv+xlf)/cp
 real(kind_phys):: xlvcp        !<= xlv/cp
 real(kind_phys):: g_inv        !<= 1./grav
! thresholds for aerosol mixing. Needed until the surface emissions are updated.
 real(kind_phys),parameter:: wfa_max = 800e10    !kg-1
 real(kind_phys),parameter:: wfa_min = 1e6       !kg-1
 real(kind_phys),parameter:: ifa_max = 500e6     !kg-1
 real(kind_phys),parameter:: ifa_min = 0.0       !kg-1
 real(kind_phys),parameter:: wfa_ht  = 2000.     !meters
 real(kind_phys),parameter:: ifa_ht  = 10000.    !meters

!-------------------------------------------------------
!---------     MYNN-specific parameters    -------------
!-------------------------------------------------------
! The parameters below depend on stability functions of module_sf_mynn.
 real(kind_phys),parameter:: cphm_st = 5.0, cphm_unst = 16.0
 real(kind_phys),parameter:: cphh_st = 5.0, cphh_unst = 16.0
! Constants for min tke in elt integration (qmin), max z/L in els (zmax),
! and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
 real(kind_phys),parameter:: qmin    = 0.0
 real(kind_phys),parameter:: zmax    = 1.0
 real(kind_phys),parameter:: Sqfac   = 3.0
 real(kind_phys),parameter:: qkemin  = 1.e-5

! Constants for cloud PDF (mym_condensation)
 real(kind_phys),parameter:: rr2     = 0.7071068
 real(kind_phys),parameter:: rrp     = 0.3989423
 
! Closure constants
 real(kind_phys),parameter:: pr      = 0.74       !prandtl number at the unstable limit
 real(kind_phys),parameter:: g1      = 0.235      !NN2009 = 0.235
 real(kind_phys),parameter:: b1      = 24.0
 real(kind_phys),parameter:: b2      = 15.0       ! CKmod    NN2009
 real(kind_phys),parameter:: c2      = 0.729      ! 0.729     0.75
 real(kind_phys),parameter:: c3      = 0.340      ! 0.340     0.352
 real(kind_phys),parameter:: c4      = 0.0
 real(kind_phys),parameter:: c5      = 0.2
 real(kind_phys),parameter:: a1      = b1*( 1.0-3.0*g1 )/6.0
 !real(kind_phys),parameter:: c1      = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) )
 real(kind_phys),parameter:: c1      = g1 -1.0/( 3.0*a1*2.88449914061481660)
 real(kind_phys),parameter:: a2      = a1*( g1-c1 )/( g1*pr )
 real(kind_phys),parameter:: g2      = b2/b1*( 1.0-c3 ) +2.0*a1/b1*( 3.0-2.0*c2 )
! parameters derived from closure constants
 real(kind_phys),parameter:: cc2     = 1.0-c2
 real(kind_phys),parameter:: cc3     = 1.0-c3
 real(kind_phys),parameter:: e1c     = 3.0*a2*b2*cc3
 real(kind_phys),parameter:: e2c     = 9.0*a1*a2*cc2
 real(kind_phys),parameter:: e3c     = 9.0*a2*a2*cc2*( 1.0-c5 )
 real(kind_phys),parameter:: e4c     = 12.0*a1*a2*cc2
 real(kind_phys),parameter:: e5c     = 6.0*a1*a1

!-------------------------------------------------------
!-------     commonly used numbers     -----------------
!-------------------------------------------------------
 real(kind_phys),parameter:: zero   = 0.0
 real(kind_phys),parameter:: one    = 1.0
 real(kind_phys),parameter:: two    = 2.0
 real(kind_phys),parameter:: three   = 3.0
 real(kind_phys),parameter:: four    = 4.0
 real(kind_phys),parameter:: five    = 5.0
 real(kind_phys),parameter:: six     = 6.0
 real(kind_phys),parameter:: seven   = 7.0
 real(kind_phys),parameter:: eight   = 8.0
 real(kind_phys),parameter:: nine    = 9.0
 real(kind_phys),parameter:: ten     = 10.0
 real(kind_phys),parameter:: twenty  = 20.0
 real(kind_phys),parameter:: thirty  = 30.0
 real(kind_phys),parameter:: forty   = 40.0
 real(kind_phys),parameter:: fifty   = 50.0
 real(kind_phys),parameter:: hundred = 100.0
 real(kind_phys),parameter:: p01     = 0.01
 real(kind_phys),parameter:: p1      = 0.1
 real(kind_phys),parameter:: p2      = 0.2
 real(kind_phys),parameter:: p25     = 0.25
 real(kind_phys),parameter:: p3      = 0.3
 real(kind_phys),parameter:: p333    = 1.0/3.0
 real(kind_phys),parameter:: p4      = 0.4
 real(kind_phys),parameter:: p5      = 0.5
 real(kind_phys),parameter:: p6      = 0.6
 real(kind_phys),parameter:: p666    = 2.0/3.0
 real(kind_phys),parameter:: p7      = 0.7
 real(kind_phys),parameter:: p75     = 0.75
 real(kind_phys),parameter:: p8      = 0.8
 real(kind_phys),parameter:: p9      = 0.9
 real(kind_phys),parameter:: p95     = 0.95
 real(kind_phys),parameter:: p99     = 0.99 

 end module module_bl_mynnedmf_common
