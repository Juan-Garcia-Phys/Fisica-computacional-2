! ============================================================
! Programa:  programa_wien
! Propósito: hallar x5 y x3 (raíces positivas de x = n(1 - e^{-x}))
!            con bisección, Newton y secante; repetir Newton desde
!            x0 = 0.5 y calcular la constante de Wien b = hc/(kB x5).
! Salida:    resultado_wien.dat, iteraciones_wien.dat, constante_wien.dat
! ============================================================
program programa_wien
  use modulo_ceros, only: dp, F, biseccion, newton, secante
  implicit none
  ! Constantes exactas del SI (2019)
  real(dp), parameter :: h_planck = 6.62607015e-34_dp    ! J s
  real(dp), parameter :: c_luz    = 299792458.0_dp       ! m/s
  real(dp), parameter :: k_boltz  = 1.380649e-23_dp      ! J/K
  real(dp), parameter :: tol = 1.0e-10_dp
  integer,  parameter :: exponentes(2) = [5, 3]          ! n = 5 (B_lambda), n = 3 (B_nu)

  real(dp) :: raiz(3, 2)          ! raiz(método, caso de n)
  integer  :: iter(3, 2)          ! iteraciones(método, caso de n)
  real(dp) :: x_trivial, b_wien
  integer  :: it_trivial, j, m, n, u_res, u_it, u_reg, u_cte
  integer  :: n_predicho

  ! Iteraciones que predice la cota de bisección (sección 2):
  ! n >= log((b-a)/eps)/log 2 con [1, 10] y eps = 1e-10
  n_predicho = ceiling(log(9.0_dp/tol)/log(2.0_dp))
  print '(a,i0)', 'Bisección: iteraciones predichas por la cota = ', n_predicho

  open(newunit=u_it, file='iteraciones_wien.dat', status='replace', action='write')

  do j = 1, 2
     n = exponentes(j)
     ! solo se registran las iteraciones del caso n = 5
     u_reg = 0
     if (n == 5) u_reg = u_it
     call biseccion(n, 1.0_dp, 10.0_dp, tol, u_reg, raiz(1, j), iter(1, j))
     call newton   (n, 10.0_dp,         tol, u_reg, raiz(2, j), iter(2, j))
     call secante  (n, 10.0_dp, 9.0_dp, tol, u_reg, raiz(3, j), iter(3, j))
  end do
  close(u_it)

  ! Newton desde x0 = 0.5: a la izquierda del punto donde F' = 0 (x = ln 5)
  call newton(5, 0.5_dp, tol, 0, x_trivial, it_trivial)

  ! Constante de Wien con x5 de Newton
  b_wien = h_planck*c_luz/(k_boltz*raiz(2, 1))

  ! --- resultado_wien.dat: n, método, raíz, iteraciones
  open(newunit=u_res, file='resultado_wien.dat', status='replace', action='write')
  do j = 1, 2
     do m = 1, 3
        write(u_res, '(i3,i3,es24.16,i5)') exponentes(j), m, raiz(m, j), iter(m, j)
     end do
  end do
  write(u_res, '(i3,i3,es24.16,i5)') 5, 4, x_trivial, it_trivial
  close(u_res)

  ! --- constante_wien.dat: x5, x3, b [m K]
  open(newunit=u_cte, file='constante_wien.dat', status='replace', action='write')
  write(u_cte, '(3es24.16)') raiz(2, 1), raiz(2, 2), b_wien
  close(u_cte)

  ! --- resumen en pantalla, con el residuo como revisión final
  print '(/,a)', '  n  método            raíz           iter    |F(raíz)|'
  do j = 1, 2
     do m = 1, 3
        print '(i3,i5,f22.14,i7,es13.3)', exponentes(j), m, raiz(m, j), iter(m, j), &
              abs(F(raiz(m, j), exponentes(j)))
     end do
  end do
  print '(i3,a,f22.14,i7,es13.3)', 5, '  N(0.5)', x_trivial, it_trivial, abs(F(x_trivial, 5))
  print '(/,a,es16.9,a)', 'Constante de Wien b = ', b_wien, ' m K'
end program programa_wien
