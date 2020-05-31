program main

use test, only: test_value
use io, only: stdout
use solver, only: solve_dsysv

use misfit

implicit none

integer :: i, j
double precision :: obs(9), pre(9), cov_diag(9), cov(9,9), rms, chi2

obs(1) = 5.79834d0
obs(2) = 28.6228d0
obs(3) = 9.57282d0
obs(4) = 39.1214d0
obs(5) = 38.0307d0
obs(6) = 21.6525d0
obs(7) = 10.7922d0
obs(8) = 43.0699d0
obs(9) = 30.9378d0

pre(1) = 3.15816d0
pre(2) = 29.0537d0
pre(3) = 9.24886d0
pre(4) = 41.3154d0
pre(5) = 39.7201d0
pre(6) = 23.4536d0
pre(7) = 9.83367d0
pre(8) = 42.3546d0
pre(9) = 33.7852d0

call misfit_rms(obs,pre,9,rms)
call test_value(rms,1.7547068607636751d0,'rms')

cov_diag(1) = 1.08d0
cov_diag(2) = 0.21d0
cov_diag(3) = 0.85d0
cov_diag(4) = 0.75d0
cov_diag(5) = 1.28d0
cov_diag(6) = 0.20d0
cov_diag(7) = 0.23d0
cov_diag(8) = 0.03d0
cov_diag(9) = 1.36d0
call misfit_chi2(obs,pre,cov_diag,.true.,9,chi2)
call test_value(chi2,59.340946552351653d0,'chi2: diagonal covariance matrix')

#ifdef USE_LAPACK
do i = 1,9
    do j = 1,9
        if (i.eq.j) then
            cov(i,j) = cov_diag(i)
        else
            cov(i,j) = (dble(mod(i,3)+mod(j,4))-3.17d0)/67.0d0
        endif
    enddo
enddo
call misfit_chi2(obs,pre,cov,.false.,9,chi2)
call test_value(chi2,61.763718745970010d0,'chi2: non-diagonal covariance matrix')
#endif

write(stdout,*) 'misfit_module unit test passed'
end
