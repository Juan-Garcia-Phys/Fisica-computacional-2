! ============================================================
! Programa:  programa_energia   (Taller 7, Parte B)
! Propósito: energía disipada por la fricción en [0, T] con
!            Romberg, E_dis = int_0^T 2 gamma v(t)^2 dt, y
!            comparación con E(0) - E(T).
! Salida:    tabla_energia.dat    (una fila por nivel: S(k,0..k))
!            resultado_energia.dat (k_final, evaluaciones, E_dis,
!                                   E(0) - E(T), diferencia)
! ============================================================
program programa_energia
  use modulo_romberg, only: dp, romberg, energia, INTEG_ENERGIA
  implicit none
  real(dp), parameter :: T   = 3.0_dp      ! tiempo final [s]
  real(dp), parameter :: tol = 1.0e-6_dp   ! tolerancia de la diagonal
  integer,  parameter :: kmax = 20         ! nivel máximo (2**20 subintervalos)
  real(dp), allocatable :: S(:,:)
  real(dp) :: e_dis, e_ref
  integer  :: k_final, n_eval, k, u
  logical  :: convergio

  call romberg(INTEG_ENERGIA, 0.0_dp, T, tol, kmax, S, k_final, n_eval, convergio)
  if (.not. convergio) print '(a)', 'Aviso: Romberg no alcanzó la tolerancia en kmax niveles'

  e_dis = S(k_final, k_final)              ! mejor valor: último elemento diagonal
  e_ref = energia(0.0_dp) - energia(T)     ! referencia exacta por conservación

  ! --- Tabla de Romberg: fila k con S(k,0), ..., S(k,k)
  open(newunit=u, file='tabla_energia.dat', status='replace', action='write')
  write(u, '(a)') '# Tabla de Romberg, E_dis en [0, 3] s. Fila k: subintervalos, S(k,0) ... S(k,k)'
  do k = 0, k_final
     write(u, '(i8, *(es22.13))') 2**k, S(k, 0:k)
  end do
  close(u)

  ! --- Resultado en una línea
  open(newunit=u, file='resultado_energia.dat', status='replace', action='write')
  write(u, '(i4, i8, 3es24.15)') k_final, n_eval, e_dis, e_ref, abs(e_dis - e_ref)
  close(u)

  ! --- Resumen en pantalla
  print '(a)', 'Tabla de Romberg (diagonal):'
  do k = 0, k_final
     print '(a,i2,a,i5,a,f16.10)', '  k =', k, '  N =', 2**k, '  S(k,k) =', S(k,k)
  end do
  print '(/,a,i3,a,i5)',  'k_final = ', k_final, '   evaluaciones = ', n_eval
  print '(a,f16.10,a)', 'E_dis (Romberg) = ', e_dis, ' J/kg'
  print '(a,f16.10,a)', 'E(0) - E(T)     = ', e_ref, ' J/kg'
  print '(a,es12.3)',   'diferencia      = ', abs(e_dis - e_ref)
end program programa_energia
