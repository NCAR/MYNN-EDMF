module module_bl_mynnedmf_diags
  use module_bl_mynnedmf_common, only: kind_phys,grav

  implicit none
  real(kind_phys),parameter::zero=0.0
  
  private
  public :: cloud_water_path !, wspd_lev, pblh, cloud_ceiling

  contains

  subroutine cloud_water_path (its, ite, kts, kte, jts, jte,                     &
                        ims, ime, kms, kme, jms, jme,                            &
                        dp, qc, qi, qs, qc_bl, qi_bl, cldfra_bl,                 &
                        lwp, iwp, swp)
    implicit none

    integer, intent(in) :: its, ite, kts, kte, jts, jte,                         &
        ims, ime, jms, jme, kms, kme
    real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in) ::           &
        dp, qc, qi, qs     
    real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in), optional :: &
        qc_bl, qi_bl, cldfra_bl     
    real(kind_phys), dimension(:, :), intent(inout), optional ::     &
        lwp, iwp, swp

    ! local variables
    real(kind_phys) :: sum1, sum2, sum3, qctotal, qitotal, qstotal
    integer :: i, k, j

    !---------------------------------------
    !Begin looping in the i- and j-direction
    !---------------------------------------

    do j = jts,jte
    do i = its,ite
       
       sum1=zero
       sum2=zero
       sum3=zero
      
      do k = kts, kte
        if (present(qc_bl) .and. present(cldfra_bl) .and. qc(i,k,j)<1e-6 .and.  & 
          cldfra_bl(i,k,j)>0.01) then
          qctotal = qc_bl(i,k,j)
        else
          qctotal = qc(i,k,j)
        endif

        if (present(qi_bl) .and. present(cldfra_bl) .and. qi(i,k,j)<1e-9 .and.  & 
          cldfra_bl(i,k,j)>0.01) then
          qitotal = qi_bl(i,k,j)
        else
          qitotal = qi(i,k,j)
        endif

        qstotal = qs(i,k,j) !there currently is no sgs snow, it's ice+snow...

        sum1 = sum1 + max((dp(i,k,j)/grav) * qctotal, zero)
        sum2 = sum2 + max((dp(i,k,j)/grav) * (qitotal+qstotal), zero) !actually frozen water path
        sum3 = sum3 + max((dp(i,k,j)/grav) * qstotal, zero)
      enddo

      lwp(i,j) = sum1 * 1000._kind_phys ! kg m-2  --> g m-2
      iwp(i,j) = sum2 * 1000._kind_phys
      swp(i,j) = sum3 * 1000._kind_phys
      print*, 'lwp', lwp

    enddo
    enddo

  end subroutine cloud_water_path

!=================================================================================================================
 end module module_bl_mynnedmf_diags
!=================================================================================================================
