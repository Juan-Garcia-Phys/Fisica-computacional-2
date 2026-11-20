! ============================================================
! Programa:  programa_sigma
! Propósito: calcular I = int_0^inf x^3/(e^x - 1) dx con el módulo
!            de Romberg del Taller 7:
!            (1) con el cambio de variable x = t/(1-t) en [0, 1] y
!            (2) sin cambio de variable en [0, 100];
!            obtener sigma = 2 pi kB^4/(h^3 c^2) I y comparar con pi^4/15.
! Salida:    resultado_sigma.dat (k_final, evaluaciones, I, sigma)
! ============================================================
program programa_sigma
  use modulo_romberg, only: dp, pi, romberg, INTEG_PLANCK_T, INTEG_PLANCK
  implicit none
  real(dp), parameter :: h_planck = 6.62607015e-34_dp    ! J s
  real(dp), parameter :: c_luz    = 299792458.0_dp       ! m/s
  real(dp), parameter :: k_boltz  = 1.380649e-23_dp      ! J/K
  real(dp), parameter :: tol = 1.0e-8_dp
  integer,  parameter :: kmax = 20                       ! límite de niveles
  real(dp), allocatable :: S(:, :)                       ! tabla de Romberg (la reserva romberg)
  real(dp) :: factor, integral(2), sigma(2), I_exacta
  integer  :: k_fin(2), n_ev(2), fila, u
  logical  :: convergio(2)

  factor   = 2.0_dp*pi*k_boltz**4/(h_planck**3*c_luz**2)
  I_exacta = pi**4/15.0_dp

  ! (1) con cambio de variable, t en [0, 1]
  call romberg(INTEG_PLANCK_T, 0.0_dp, 1.0_dp, tol, kmax, S, k_fin(1), n_ev(1), convergio(1))
  integral(1) = S(k_fin(1), k_fin(1))
  ! (2) sin cambio de variable, x en [0, 100]
  !     (S es intent(out) allocatable: se libera y se reserva de nuevo)
  call romberg(INTEG_PLANCK, 0.0_dp, 100.0_dp, tol, kmax, S, k_fin(2), n_ev(2), convergio(2))
  integral(2) = S(k_fin(2), k_fin(2))
  sigma = factor*integral

  do fila = 1, 2
     if (.not. convergio(fila)) print '(a,i0,a)', 'fila ', fila, ': Romberg no convergió en kmax niveles'
  end do

  open(newunit=u, file='resultado_sigma.dat', status='replace', action='write')
  do fila = 1, 2
     write(u, '(i4,i10,2es24.16)') k_fin(fila), n_ev(fila), integral(fila), sigma(fila)
  end do
  close(u)

  print '(a)', ' fila  k_final  evaluaciones        I                 error |I - pi^4/15|'
  do fila = 1, 2
     print '(i4,i9,i12,es24.15,es14.3)', fila, k_fin(fila), n_ev(fila), integral(fila), &
           abs(integral(fila) - I_exacta)
  end do
  print '(/,a,es24.15)', 'pi^4/15             = ', I_exacta
  print '(a,es16.9,a)',  'sigma (con cambio)  = ', sigma(1), ' W m^-2 K^-4'
end program programa_sigma
