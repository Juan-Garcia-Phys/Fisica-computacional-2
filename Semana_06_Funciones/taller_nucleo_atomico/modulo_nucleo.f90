! Taller 6 -- Módulo de la fórmula semiempírica de masa (modelo de la gota líquida)
! Contiene las constantes del modelo y todos los procedimientos de cálculo.
! El programa principal (programa_nucleo.f90) lo usa con "use modulo_nucleo, only: ...".
! Convención: c = 1, energías y masas en MeV; E_B < 0 para un núcleo ligado.
module modulo_nucleo
  implicit none

  ! --- Constantes globales del modelo (todas parameter: R3) ---
  integer,  parameter :: dp  = kind(1.0d0)     ! reales de doble precisión (R7)
  real(dp), parameter :: a_v = 15.8_dp         ! coeficiente de volumen [MeV]
  real(dp), parameter :: a_s = 18.3_dp         ! coeficiente de superficie [MeV]
  real(dp), parameter :: a_c = 0.71_dp         ! coeficiente de Coulomb [MeV]
  real(dp), parameter :: a_a = 92.7_dp         ! coeficiente de asimetría [MeV]
  real(dp), parameter :: m_n = 939.57_dp       ! masa del neutrón [MeV]
  real(dp), parameter :: m_h = 938.78_dp       ! masa del átomo de hidrógeno [MeV]

contains

  ! ------------------------------------------------------------------
  ! e_volumen: término de volumen de la energía de enlace por nucleón
  !   Entradas: ninguna (por nucleón es constante: no depende de Z ni de A)
  !   Salida:   -a_V [MeV]
  ! ------------------------------------------------------------------
  real(dp) function e_volumen()
    e_volumen = -a_v
  end function e_volumen

  ! ------------------------------------------------------------------
  ! e_superficie: término de superficie por nucleón, a_S * A^(-1/3)
  !   Entrada: a (número de nucleones, entero > 0)
  !   Salida:  término de superficie [MeV]
  ! ------------------------------------------------------------------
  real(dp) function e_superficie(a)
    integer, intent(in) :: a
    ! real(a, dp) y el exponente real evitan la división entera -1/3 = 0
    e_superficie = a_s*real(a, dp)**(-1.0_dp/3.0_dp)
  end function e_superficie

  ! ------------------------------------------------------------------
  ! e_coulomb: término de Coulomb por nucleón, a_C * Z^2 * A^(-4/3)
  !   Entradas: z (número de protones), a (número de nucleones)
  !   Salida:   término de Coulomb [MeV]
  ! ------------------------------------------------------------------
  real(dp) function e_coulomb(z, a)
    integer, intent(in) :: z, a
    e_coulomb = a_c*real(z, dp)**2*real(a, dp)**(-4.0_dp/3.0_dp)
  end function e_coulomb

  ! ------------------------------------------------------------------
  ! e_asimetria: término de asimetría por nucleón, a_A * (Z - A/2)^2 / A^2
  !   Entradas: z (número de protones), a (número de nucleones)
  !   Salida:   término de asimetría [MeV]
  ! ------------------------------------------------------------------
  real(dp) function e_asimetria(z, a)
    integer, intent(in) :: z, a
    real(dp) :: exceso                           ! local: Z - A/2, calculado en reales
    ! A/2 se calcula en reales: con enteros, A impar perdería el 0.5
    exceso      = real(z, dp) - real(a, dp)/2.0_dp
    e_asimetria = a_a*exceso**2/real(a, dp)**2
  end function e_asimetria

  ! ------------------------------------------------------------------
  ! energia_enlace: energía de enlace por nucleón E_B(Z, A)
  !   Entradas: z (protones), a (nucleones)
  !   Salida:   E_B [MeV], suma de los cuatro términos (negativa si está ligado)
  ! ------------------------------------------------------------------
  real(dp) function energia_enlace(z, a)
    integer, intent(in) :: z, a
    ! Reutiliza las cuatro funciones: la fórmula de cada término está en un solo lugar
    energia_enlace = e_volumen() + e_superficie(a) + e_coulomb(z, a) + e_asimetria(z, a)
  end function energia_enlace

  ! ------------------------------------------------------------------
  ! masa_atomica: masa atómica aproximada M(Z, A) = Z*M_H + (A-Z)*M_n + A*E_B
  !   Entradas: z (protones), a (nucleones)
  !   Salida:   M [MeV] (energía equivalente, c = 1)
  ! ------------------------------------------------------------------
  real(dp) function masa_atomica(z, a)
    integer, intent(in) :: z, a
    ! A*E_B es la energía de enlace total: negativa, reduce la masa
    masa_atomica = real(z, dp)*m_h + real(a - z, dp)*m_n + real(a, dp)*energia_enlace(z, a)
  end function masa_atomica

  ! ------------------------------------------------------------------
  ! z_estable: para un A dado, busca el Z de máximo enlace (menor E_B)
  !   Entrada:  a (número de nucleones)
  !   Salidas:  z_min (Z con el menor E_B entre Z = 1, ..., A)
  !             eb_min (ese E_B mínimo [MeV])
  ! Es una subrutina porque entrega DOS resultados separados (regla 3.3):
  ! una función solo devuelve un valor por su nombre.
  ! ------------------------------------------------------------------
  subroutine z_estable(a, z_min, eb_min)
    integer,  intent(in)  :: a
    integer,  intent(out) :: z_min
    real(dp), intent(out) :: eb_min
    integer  :: z                                ! local: candidato a número de protones
    real(dp) :: eb                               ! local: E_B del candidato

    ! El primer candidato (Z = 1) es el mínimo provisional; se asigna en el
    ! cuerpo (no en la declaración) para que cada llamada empiece de cero.
    z_min  = 1
    eb_min = energia_enlace(1, a)
    do z = 2, a
       eb = energia_enlace(z, a)
       if (eb < eb_min) then                     ! más negativo = más ligado
          z_min  = z
          eb_min = eb
       end if
    end do
  end subroutine z_estable

end module modulo_nucleo
