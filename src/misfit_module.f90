!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------- MISFIT ROUTINES ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
module misfit

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine misfit_rms(obs,pre,n,rms)
!----
! Calculate the root-mean-squared misfit between an observed and predicted set of points:
!
!     RMS = sqrt(sum((obs-pre)^2)/n)
!-----

use io, only: stderr

implicit none

! Arguments
integer :: n
double precision :: obs(n), pre(n), rms

! Local variables
integer :: i
double precision :: dif(n)


! Check that n is greater than zero
if (n.le.0) then
    write(stderr,*) 'misfit_rms: n is less than or equal to zero'
endif

! Initialize RMS misfit
rms = 0.0d0

! Add differences squared to RMS
do i = 1,n
    dif(i) = obs(i)-pre(i)
    dif(i) = dif(i)*dif(i)
    rms = rms + dif(i)
enddo

! Take square root
rms = sqrt(rms/n)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine misfit_chi2(obs,pre,cov,isCovMatrixDiagonal,n,chi2)
!----
! Calculate the chi-squared misfit between an observed and predicted set of points:
!
!     chi2 = dif_tr * cov_mat^(-1) * dif
!
! where
!
!     dif = obs-pre
!     dif_tr = transpose(dif)
!     cov_mat = covariance matrix
!-----

#ifdef USE_LAPACK
use solver, only: solve_dsysv
#endif

use io, only: stderr
use random, only: iseed, r8_normal_ab

implicit none

! Arguments
integer :: n
double precision :: obs(n), pre(n), cov(n,n), chi2, model_uncertainty
logical :: isCovMatrixDiagonal

! Local variables
integer :: i, ierr
double precision :: dif(n), vec(n)


! Check that n is greater than zero
if (n.le.0) then
    write(stderr,*) 'misfit_chi2: n is less than or equal to zero'
endif

! Initialize chi-squared value
chi2 = 0.0d0

! Calculate difference between observed and predicted
dif = obs-pre

! Calculate dif_trans*cov^(-1)*dif
if (isCovMatrixDiagonal) then
    do i = 1,n
        ! if (model_uncertainty.le.0.0d0) then
            chi2 = chi2 + dif(i)*dif(i)/cov(i,1)
        ! else
            ! dif(i) = dif(i) - sign(model_uncertainty*obs(i),dif(i))
            ! dif(i) = dif(i) + r8_normal_ab(0.0d0,model_uncertainty*obs(i),iseed)
            ! chi2 = chi2 + dif(i)*dif(i)/ (cov(i,1) + model_uncertainty**2 * obs(i)**2)
            ! chi2 = chi2 + dif(i)*dif(i)/ (cov(i,1))
        ! endif
    enddo
else

#ifdef USE_LAPACK
    call solve_dsysv(cov,dif,vec,n,ierr)
    do i = 1,n
        chi2 = chi2 + dif(i)*vec(i)
    enddo
#else
    call usage('misfit_chi2: cannot calculate chi-squared for non-diagonal covariance matrix '//&
               'without LAPACK')
#endif

endif

return
end subroutine

end module
