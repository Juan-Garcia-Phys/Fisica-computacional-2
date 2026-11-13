! ============================================================
! Programa:  programa_pi   (Taller 7, Parte D)
! Propósito: pi = int_0^1 4/(1+x^2) dx con el mismo módulo de
!            Romberg de la Parte B, sin modificarlo.
! Salida:    resultado_pi.dat (k_final, evaluaciones, valor, error)
! ============================================================
program programa_pi
  use modulo_romberg, only: dp, pi, romberg, INTEG_PI
  implicit none
  real(dp), parameter :: tol = 1.0e-5_dp
  integer,  parameter :: kmax = 20
  real(dp), allocatable :: S(:,:)
  real(dp) :: valor
  integer  :: k_final, n_eval, k, u
  logical  :: convergio

  call romberg(INTEG_PI, 0.0_dp, 1.0_dp, tol, kmax, S, k_final, n_eval, convergio)
  if (.not. convergio) print '(a)', 'Aviso: Romberg no alcanzó la tolerancia en kmax niveles'
  valor = S(k_final, k_final)

  open(newunit=u, file='resultado_pi.dat', status='replace', action='write')
  write(u, '(i4, i8, 2es24.15)') k_final, n_eval, valor, abs(valor - pi)
  close(u)

  print '(a)', 'Diagonal de la tabla de Romberg para pi:'
  do k = 0, k_final
     print '(a,i2,a,f14.10,a,es10.2)', '  k =', k, '  S(k,k) =', S(k,k), '  error =', abs(S(k,k) - pi)
  end do
  print '(/,a,i3,a,i5)', 'k_final = ', k_final, '   evaluaciones = ', n_eval
  print '(a,f14.10,a,es10.2)', 'pi (Romberg) = ', valor, '   error = ', abs(valor - pi)
end program programa_pi
