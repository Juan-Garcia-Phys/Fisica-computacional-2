! ============================================================
! Módulo:    modulo_romberg
! Propósito: integración de Romberg (tabla triangular con
!            extrapolación de Richardson) con conteo de
!            evaluaciones, y los integrandos del Taller 7:
!              opción 1 -> 2*gamma*v(t)**2   (energía disipada, Parte B)
!              opción 2 -> 4/(1 + x**2)      (valor de pi, Parte D)
!            Se usa sin cambios en las Partes B y D.
! ============================================================
module modulo_romberg
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 4*atan(1.0_dp)

  ! Parámetros del oscilador amortiguado x(t) = exp(-gamma t) cos(omega t)
  real(dp), parameter :: gamma = 0.5_dp           ! amortiguamiento [1/s]
  real(dp), parameter :: omega = 2*pi             ! frecuencia angular [rad/s]
  real(dp), parameter :: omega0_2 = omega**2 + gamma**2   ! omega_0^2 [1/s^2]

  ! Selección del integrando
  integer, parameter :: INTEG_ENERGIA = 1
  integer, parameter :: INTEG_PI      = 2

contains

  ! ----------------------------------------------------------
  ! posicion: x(t) del oscilador amortiguado [m]
  ! ----------------------------------------------------------
  elemental real(dp) function posicion(t)
    real(dp), intent(in) :: t
    posicion = exp(-gamma*t)*cos(omega*t)
  end function posicion

  ! ----------------------------------------------------------
  ! velocidad: v(t) exacta del oscilador [m/s]
  ! ----------------------------------------------------------
  elemental real(dp) function velocidad(t)
    real(dp), intent(in) :: t
    velocidad = -exp(-gamma*t)*(gamma*cos(omega*t) + omega*sin(omega*t))
  end function velocidad

  ! ----------------------------------------------------------
  ! energia: E(t) = v^2/2 + omega0^2 x^2/2, por unidad de masa [J/kg]
  ! ----------------------------------------------------------
  elemental real(dp) function energia(t)
    real(dp), intent(in) :: t
    energia = 0.5_dp*velocidad(t)**2 + 0.5_dp*omega0_2*posicion(t)**2
  end function energia

  ! ----------------------------------------------------------
  ! integrando: integrando elegido por 'opcion' (elemental:
  ! acepta un arreglo de puntos)
  ! Entradas: opcion (INTEG_ENERGIA o INTEG_PI), x
  ! Salida:   f(x); con una opción desconocida devuelve 0
  ! ----------------------------------------------------------
  elemental real(dp) function integrando(opcion, x)
    integer,  intent(in) :: opcion
    real(dp), intent(in) :: x
    select case (opcion)
    case (INTEG_ENERGIA)
       integrando = 2*gamma*velocidad(x)**2      ! potencia disipada dE/dt con signo +
    case (INTEG_PI)
       integrando = 4/(1 + x**2)
    case default
       integrando = 0.0_dp
    end select
  end function integrando

  ! ----------------------------------------------------------
  ! romberg: integral de 'opcion' en [a, b] por la tabla de Romberg
  ! Entradas: opcion, a, b, tol (tolerancia del criterio de la
  !           diagonal), kmax (nivel máximo)
  ! Salidas:  S(0:kmax, 0:kmax), tabla triangular (S(k,j) con j<=k);
  !           k_final, último nivel calculado;
  !           n_eval, evaluaciones del integrando (2**k_final + 1);
  !           convergio, .true. si |S(k,k) - S(k-1,k-1)| < tol
  ! ----------------------------------------------------------
  subroutine romberg(opcion, a, b, tol, kmax, S, k_final, n_eval, convergio)
    integer,  intent(in)  :: opcion, kmax
    real(dp), intent(in)  :: a, b, tol
    real(dp), allocatable, intent(out) :: S(:,:)
    integer,  intent(out) :: k_final, n_eval
    logical,  intent(out) :: convergio
    real(dp) :: h, suma, factor
    integer  :: k, j, i, m

    ! Tabla con índices desde 0, como en las fórmulas
    allocate(S(0:kmax, 0:kmax))
    S = 0.0_dp

    ! Nivel 0: un solo trapecio sobre todo el intervalo
    h = b - a
    S(0,0) = 0.5_dp*h*(integrando(opcion, a) + integrando(opcion, b))
    n_eval = 2
    convergio = .false.
    k_final = kmax

    do k = 1, kmax
       h = 0.5_dp*h
       ! Solo los puntos nuevos (puntos medios del nivel anterior):
       ! x = a + (2i - 1) h, i = 1, ..., 2**(k-1)
       m = 2**(k - 1)
       suma = sum(integrando(opcion, [( a + (2*i - 1)*h, i = 1, m )]))
       n_eval = n_eval + m
       S(k,0) = 0.5_dp*S(k-1,0) + h*suma

       ! Extrapolación de Richardson (4**j en real: evita desborde)
       do j = 1, k
          factor = 4.0_dp**j
          S(k,j) = (factor*S(k,j-1) - S(k-1,j-1))/(factor - 1)
       end do

       ! Criterio de parada: dos elementos diagonales consecutivos
       if (abs(S(k,k) - S(k-1,k-1)) < tol) then
          k_final = k
          convergio = .true.
          exit
       end if
    end do
  end subroutine romberg

end module modulo_romberg
