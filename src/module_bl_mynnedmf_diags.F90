module module_bl_mynnedmf_diags
  use module_bl_mynnedmf_common, only: kind_phys,g_inv

  implicit none
  real(kind_phys),parameter::zero=0.0
  
!  private
!  public :: cloud_water_path, wspd_at_hgts, cloud_ceiling ! pblh

  contains

!===================================================================
! Subroutine to call MYNN-EDMF diagnostics
!===================================================================
  subroutine mynnedmf_diags (&
                kts     , kte    , delp1     , dz1 , zw1   , zagl1, u1  , v1  , tk1 , qc1 , qi1 , qs1,&
                qc_bl1  , qi_bl1 , cldfra_bl1, rho1, xland1, pblh1, lwp1, iwp1, swp1, cldceil1       ,&
                wspd101 , wspd801, wspd1601  , maxcldfra1  , maxcldfra_pbl1   , bl_mynn_diags  )

    implicit none

    integer, intent(in) :: kts, kte
    integer, intent(in) :: bl_mynn_diags
    real(kind_phys), intent(in) :: xland1, pblh1
    real(kind_phys), dimension(kts:kte), intent(in) :: delp1, dz1, zagl1, u1, v1, tk1, qc1, &
                       qi1, qs1, rho1
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw1
    real(kind_phys), dimension(kts:kte), intent(in), optional :: qc_bl1, qi_bl1,     &
                       cldfra_bl1 
    real(kind_phys), dimension(kts:kte) :: qctotal1, qitotal1, qstotal1
    real(kind_phys), intent(inout), optional :: lwp1, iwp1, swp1, maxcldfra1, maxcldfra_pbl1
    real(kind_phys), intent(out)  , optional :: cldceil1, wspd101, wspd801, wspd1601

    integer :: k

    do k = kts, kte
      if (qc1(k)<1e-6 .and. cldfra_bl1(k)>0.01) then
        qctotal1(k) = qc_bl1(k)
      else
        qctotal1(k) = qc1(k)
      endif

      if (qi1(k)<1e-9 .and. cldfra_bl1(k)>0.01) then
        qitotal1(k) = qi_bl1(k)
      else
        qitotal1(k) = qi1(k)
      endif

      qstotal1(k) = qs1(k) !there currently is no sgs snow, it's ice+snow...
    end do

    call cloud_water_path (kts, kte, delp1, qctotal1, qitotal1, qstotal1,             &
          lwp1, iwp1, swp1)

    call cloud_ceiling (kts, kte, dz1, zagl1, zw1, tk1, qctotal1, qitotal1, qstotal1, &
          cldfra_bl1, rho1, xland1, cldceil1)

    if (bl_mynn_diags >= 2) then
      call wspd_at_hgts (kts, kte, zagl1, u1, v1, wspd101, wspd801, wspd1601)

      maxcldfra1     = maxval(cldfra_bl1)
      maxcldfra_pbl1 = maxval(cldfra_bl1, mask=(zagl1 < (pblh1 + 1000._kind_phys)))
    endif

  end subroutine mynnedmf_diags

!===================================================================
! Subroutine to calculate LWP, IWP, and SWP
!===================================================================
  subroutine cloud_water_path (kts, kte, delp1, qctotal1, qitotal1, qstotal1,        &
                lwp1, iwp1, swp1)
    implicit none

    integer, intent(in) :: kts, kte
    real(kind_phys), dimension(kts:kte), intent(in) :: delp1, qctotal1, qitotal1, qstotal1
    real(kind_phys), intent(inout), optional :: lwp1, iwp1, swp1

    ! local variables
    real(kind_phys) :: sum1, sum2, sum3
    integer :: k

    sum1=zero
    sum2=zero
    sum3=zero

    do k = kts, kte
      sum1 = sum1 + max((delp1(k) * g_inv) * qctotal1(k), zero)
      sum2 = sum2 + max((delp1(k) * g_inv) * qitotal1(k), zero)
      sum3 = sum3 + max((delp1(k) * g_inv) * qstotal1(k), zero)
    enddo

    lwp1 = sum1 * 1000._kind_phys ! kg m-2  --> g m-2
    iwp1 = sum2 * 1000._kind_phys
    swp1 = sum3 * 1000._kind_phys
    ! print*, 'lwp', lwp1

  end subroutine cloud_water_path

!===================================================================
! Subroutine to calculate wind speed at specific heights
!===================================================================
  subroutine wspd_at_hgts (kts, kte, zagl1, u1, v1, wspd101, wspd801, wspd1601)
    implicit none

    integer, intent(in) :: kts, kte
    real(kind_phys), dimension(kts:kte), intent(in) :: zagl1, u1, v1

    real(kind_phys) :: depth

    real(kind_phys), intent(out) :: wspd101, wspd801, wspd1601
    integer :: k

    wspd101   = -99._kind_phys
    wspd801   = -99._kind_phys
    wspd1601  = -99._kind_phys
    do k = kts, kte-1
      depth = zagl1(k+1) - zagl1(k)
      !print*, 'depth', depth,'zagl',zagl(k), 'wspd', sqrt(u1(k)**2+v1(k)**2)

      ! 10-m wind speed
      if (zagl1(k) > 10. .and. wspd101 == -99.) then
        wspd101 = interpolate_wind (k, kts, 10.0, zagl1(k), u1(k), v1(k),            &
                    u1(k-1), v1(k-1), depth)
        ! print*, 'k', k, 'wspd10', wspd101
      endif

      ! 80-m wind speed
      if (zagl1(k) > 80. .and. wspd801 == -99.) then
        wspd801 = interpolate_wind (k, kts, 80.0, zagl1(k), u1(k), v1(k),            &
                    u1(k-1), v1(k-1), depth)
        ! print*, 'k', k, 'wspd80', wspd801
      endif

      ! 160-m wind speed
      if (zagl1(k) > 160. .and. wspd1601 == -99.) then
        wspd1601 = interpolate_wind (k, kts, 160.0, zagl1(k), u1(k), v1(k),          &
                    u1(k-1), v1(k-1), depth)
        ! print*, 'k', k, 'wspd160', wspd1601
      endif

    enddo

  end subroutine wspd_at_hgts

!===================================================================
! Subroutine to calculate cloud ceiling based on cloud optical depth
!===================================================================
  subroutine cloud_ceiling (kts, kte, dz1, zagl1, zw1, tk1, qctotal1, qitotal1, qstotal1,     &
                cldfra_bl1, rho1, xland1, cldceil1)
    implicit none

    integer, intent(in) :: kts, kte
    real(kind_phys), intent(in) :: xland1
    real(kind_phys), dimension(kts:kte), intent(in) :: dz1, zagl1, zw1, tk1, rho1, cldfra_bl1
    real(kind_phys), dimension(kts:kte), intent(in) :: qctotal1, qitotal1, qstotal1
    real(kind_phys) :: cld_depth, cldtau1, re_i, re_c, sum_cldtau
    logical :: check_fog_layer
    real(kind_phys), intent(out) :: cldceil1
    integer :: k

    !---constants---
    real(kind_phys),parameter :: cldtau_thld  = 2.0           ! threshold for cloud optical depth
    real(kind_phys),parameter :: min_cloud_dz = 40.0          ! threshold for meaningful fog layer
    real(kind_phys),parameter :: rho_w        = 1000.0        ! water density
    real(kind_phys),parameter :: rho_i        = 917.0         ! ice density

    ! use constant water cloud effective radius over land and water
    if (xland1 == 1.0) then     ! over land
      re_c = 7.0e-6_kind_phys
    else                        ! over water
      re_c = 1.1e-5_kind_phys
    endif

    sum_cldtau = zero
    check_fog_layer = .false.
    cldceil1   = -99._kind_phys ! default value for no cloud ceiling
    cld_depth  = zero

    do k = kts, kte-1
      ! temperature-dependent ice cloud effective radius (Mishara et al. 2014)
      re_i = 173.46_kind_phys + 2.14_kind_phys * (tk1(k) - 273.15_kind_phys)

      cldtau1= (1.5_kind_phys * qctotal1(k) * rho1(k) * dz1(k) / (rho_w * re_c)) +   &
                (1.5_kind_phys * (qitotal1(k)+qstotal1(k)) * rho1(k) * dz1(k) / (rho_i *re_i))
      !print*, 'cldtau1', cldtau1, 'cldfra_bl1', cldfra_bl1(k), 'qctotal', qctotal

      if (cldtau1 > zero .and. cldfra_bl1(k) >= 0.5_kind_phys) then
        cld_depth =cld_depth + dz1(k)                             ! track cloud depth
        sum_cldtau = sum_cldtau + cldtau1

        if (k==1 .or. k==2) check_fog_layer = .true.
        if (cld_depth >= min_cloud_dz) check_fog_layer = .false.

      else
        if (check_fog_layer .and. cld_depth < min_cloud_dz) then  ! fog layer detected
          sum_cldtau = zero
          cld_depth  = zero
          check_fog_layer = .false.
          cycle ! to next k level
        endif

        ! clear sky or cloud ended; reset cloud depth for higher clouds
        sum_cldtau      = zero
        cld_depth       = zero
        check_fog_layer = .false.
      endif

      if (.not. check_fog_layer) then
        if (sum_cldtau >= cldtau_thld .and. cldfra_bl1(k) >= 0.5) then
          cldceil1 = zw1(k) + ((sum_cldtau-cldtau_thld)/cldtau1) * dz1(k)
          ! print*, 'k', k, 'cldceil1', cldceil1
          exit
        endif
      endif

    enddo

  end subroutine cloud_ceiling

  function interpolate_wind(k_idx, kts, target_z, zagl,                      &
                            u1_curr,v1_curr, u1_prev, v1_prev, depth)         &
                              result(wspd1_hgt)

    integer, intent(in) :: k_idx, kts
    real(kind_phys), intent(in)    :: target_z, zagl, u1_curr, v1_curr,      &
            u1_prev, v1_prev, depth
    real :: wspd1_hgt

    real :: wgt, wspd1_curr, wspd1_prev

    wspd1_curr = sqrt(u1_curr**2 + v1_curr**2)

    if (k_idx == kts) then
      wspd1_hgt = wspd1_curr
    else
      wgt = (zagl - target_z) /depth
      wgt = min(max(zero, wgt), 1.0_kind_phys)
      wspd1_prev = sqrt(u1_prev**2 + v1_prev**2)
      wspd1_hgt = (1.0_kind_phys - wgt) * wspd1_curr + wgt * wspd1_prev
    endif

  end function interpolate_wind

!=================================================================================================================
 end module module_bl_mynnedmf_diags
!=================================================================================================================
