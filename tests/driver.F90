! ###########################################################################################
! CI test driver for MYNN-EDMF scheme
! ###########################################################################################
program driver

use module_bl_mynnedmf_driver
use module_bl_mynnedmf_tests

implicit none

call mynnedmf_test(case='clr',bl_mynn_closure=2.6,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2,          &
      bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1,   &
      bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,&
      tke_budget=1,bl_mynn_diags=0)

call mynnedmf_test(case='clr',bl_mynn_closure=3.0,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2,          &
      bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1,   &
      bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,&
      tke_budget=1,bl_mynn_diags=0,mix_chem=.true.,enh_mix=.true.)

call mynnedmf_test(case='marine_StCu',bl_mynn_closure=3.0,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2,  &
      bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1,   &
      bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,&
      tke_budget=1,bl_mynn_diags=2)

! restart test
call restart_test(case='clr',bl_mynn_closure=2.6,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2, &
     bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1, &
     bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,&
     tke_budget=1,bl_mynn_diags=0,n_restart_in=10,tol_in=1.0e-6)

end program driver
