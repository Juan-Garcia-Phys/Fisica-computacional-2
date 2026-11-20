! ============================================================
! Módulo:    modulo_ceros
! Propósito: ecuación de Wien F(x) = x - n(1 - e^{-x}), su derivada
!            y tres métodos para hallar ceros: bisección, Newton y
!            secante. Cada método cuenta sus iteraciones y, si se le
!            pasa una unidad de archivo abierta (u_reg /= 0), registra
!            cada valor nuevo para analizar la convergencia. Se compara
!            con /= 0 y no con > 0 porque open(newunit=...) entrega
!            números de unidad NEGATIVOS.
! Convención del registro: una fila por iteración con
!            (método, número de iteración, x), numerada desde 1.
! ============================================================
module modulo_ceros
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: max_iter = 100      ! requisito 3 del taller

contains

  ! F: función cuyo cero se busca. n es el exponente del espectro
  !    (5 para B_lambda, 3 para B_nu); entra con intent(in) para que
  !    ningún método pueda modificarlo por accidente.
  real(dp) function F(x, n)
    real(dp), intent(in) :: x
    integer,  intent(in) :: n
    F = x - real(n, dp)*(1.0_dp - exp(-x))
  end function F

  ! dF: derivada analítica, F'(x) = 1 - n e^{-x}. Se anula en x = ln(n),
  !     punto que Newton debe evitar (denominador nulo).
  real(dp) function dF(x, n)
    real(dp), intent(in) :: x
    integer,  intent(in) :: n
    dF = 1.0_dp - real(n, dp)*exp(-x)
  end function dF

  ! ----------------------------------------------------------
  ! Bisección (sección 2 del notebook).
  ! Entrada: n, intervalo [a0, b0], tolerancia tol, unidad u_reg
  !          (0 = no registrar).
  ! Salida:  raiz = último punto medio, iter = iteraciones usadas.
  ! a0 y b0 se copian a variables locales para no alterar los datos
  ! iniciales del programa principal.
  ! ----------------------------------------------------------
  subroutine biseccion(n, a0, b0, tol, u_reg, raiz, iter)
    integer,  intent(in)  :: n, u_reg
    real(dp), intent(in)  :: a0, b0, tol
    real(dp), intent(out) :: raiz
    integer,  intent(out) :: iter
    real(dp) :: a, b, c
    integer  :: k

    a = a0
    b = b0
    ! Paso 1: sin cambio de signo no hay garantía de raíz en [a, b]
    if (F(a, n)*F(b, n) >= 0.0_dp) stop 'biseccion: F(a) y F(b) no tienen signos distintos'

    do k = 1, max_iter
       c = 0.5_dp*(a + b)                          ! paso 2: punto medio
       if (u_reg /= 0) write(u_reg, '(i3,i5,es24.16)') 1, k, c
       if (b - c <= tol) exit                      ! paso 3: error garantizado <= tol
       if (F(a, n)*F(c, n) <= 0.0_dp) then         ! paso 4: conservar la mitad
          b = c                                    !         con cambio de signo
       else
          a = c
       end if
    end do
    if (k > max_iter) then                         ! el do terminó sin exit
       print '(a)', 'biseccion: no convergió en max_iter iteraciones'
       k = max_iter
    end if
    raiz = c
    iter = k
  end subroutine biseccion

  ! ----------------------------------------------------------
  ! Newton-Raphson (sección 3): x_{k+1} = x_k - F(x_k)/F'(x_k).
  ! Parada: |x_{k+1} - x_k| < tol, que estima el error de x_k.
  ! ----------------------------------------------------------
  subroutine newton(n, x0, tol, u_reg, raiz, iter)
    integer,  intent(in)  :: n, u_reg
    real(dp), intent(in)  :: x0, tol
    real(dp), intent(out) :: raiz
    integer,  intent(out) :: iter
    real(dp) :: x, paso, pendiente
    integer  :: k

    x = x0
    do k = 1, max_iter
       pendiente = dF(x, n)
       if (pendiente == 0.0_dp) stop 'newton: derivada nula, la tangente no corta el eje'
       paso = -F(x, n)/pendiente                   ! corrección de Newton
       x = x + paso
       if (u_reg /= 0) write(u_reg, '(i3,i5,es24.16)') 2, k, x
       if (abs(paso) < tol) exit
    end do
    if (k > max_iter) then
       print '(a)', 'newton: no convergió en max_iter iteraciones'
       k = max_iter
    end if
    raiz = x
    iter = k
  end subroutine newton

  ! ----------------------------------------------------------
  ! Secante (sección 4.1): la derivada se reemplaza por la pendiente
  ! de la recta que pasa por las dos últimas estimaciones.
  ! La iteración k produce x_{k+1}; la primera produce x_2.
  ! Se guarda F(x_anterior) para hacer una sola evaluación nueva por paso.
  ! ----------------------------------------------------------
  subroutine secante(n, x0, x1, tol, u_reg, raiz, iter)
    integer,  intent(in)  :: n, u_reg
    real(dp), intent(in)  :: x0, x1, tol
    real(dp), intent(out) :: raiz
    integer,  intent(out) :: iter
    real(dp) :: x_ant, x_act, x_nuevo, f_ant, f_act, denominador
    integer  :: k

    x_ant = x0
    x_act = x1
    f_ant = F(x_ant, n)
    f_act = F(x_act, n)
    do k = 1, max_iter
       denominador = f_act - f_ant
       if (denominador == 0.0_dp) stop 'secante: F(x_k) = F(x_{k-1}), secante horizontal'
       x_nuevo = x_act - f_act*(x_act - x_ant)/denominador
       if (u_reg /= 0) write(u_reg, '(i3,i5,es24.16)') 3, k, x_nuevo
       ! desplazar la ventana de dos puntos antes de la prueba de parada
       x_ant = x_act
       f_ant = f_act
       x_act = x_nuevo
       f_act = F(x_act, n)
       if (abs(x_act - x_ant) < tol) exit
    end do
    if (k > max_iter) then
       print '(a)', 'secante: no convergió en max_iter iteraciones'
       k = max_iter
    end if
    raiz = x_act
    iter = k
  end subroutine secante

end module modulo_ceros
