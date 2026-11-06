! Taller 6 -- Cálculos adicionales de la sección E (preguntas 2, 3 y 4)
! Reutiliza las funciones del módulo; no modifica el modelo ni los archivos de D.
program calculos_E
  use modulo_nucleo, only: dp, energia_enlace, masa_atomica, z_estable
  implicit none
  real(dp), parameter :: a_p = 12.0_dp     ! coeficiente de apareamiento [MeV] (solo para E3)
  integer  :: z56                          ! Z de máximo enlace para A = 56
  real(dp) :: eb56, eb_fe, et_mod, et_fe, dp_mod, dp_fe
  real(dp) :: eb_max, eb_62, eb_ni, q

  ! --- E2: máximo del modelo frente al punto A = 62 de la curva y al Ni-62 ---
  eb_max = energia_enlace(29, 65)          ! máximo de la curva (maximo.dat: A=65, Z=29)
  call z_estable(62, z56, eb_62)           ! punto de la curva con A = 62 (reutiliza z56 como salida)
  eb_ni = energia_enlace(28, 62)           ! modelo evaluado en el núclido Ni-62
  print '(a)', '--- E2 ---'
  print '(a, f14.8)', 'E_B maximo del modelo (A=65, Z=29)       =', eb_max
  print '(a, i3, a, f14.8)', 'E_B curva en A=62 (Z=', z56, ')              =', eb_62
  print '(a, f14.8)', 'Diferencia E_B(max) - E_B(curva, A=62)    =', eb_max - eb_62
  print '(a, f14.6)', 'Porcentaje respecto de |E_B(curva, A=62)| =', 100.0_dp*abs(eb_max - eb_62)/abs(eb_62)
  print '(a, f14.8)', 'Modelo en Ni-62 (Z=28, A=62)              =', eb_ni
  print '(a, f14.8)', 'B/A modelo Ni-62 - B/A experimental        =', -eb_ni - 8.7945_dp

  ! --- E3: A = 56, Z del modelo frente al hierro (Z = 26) ---
  call z_estable(56, z56, eb56)
  eb_fe  = energia_enlace(26, 56)
  et_mod = 56.0_dp*eb56                    ! energía total de enlace [MeV]
  et_fe  = 56.0_dp*eb_fe
  dp_mod = apareamiento(z56, 56)
  dp_fe  = apareamiento(26, 56)
  print '(a)', '--- E3 ---'
  print '(a, i3, a, f14.8, a, f14.6)', 'Modelo: Z =', z56, '  E_B =', eb56, '  E_total =', et_mod
  print '(a, i3, a, f14.8, a, f14.6)', 'Hierro: Z =', 26,  '  E_B =', eb_fe, '  E_total =', et_fe
  print '(a, f14.6)', 'E_total(Z_modelo) - E_total(Fe)           =', et_mod - et_fe
  print '(a, f10.6, a, f10.6)', 'Apareamiento total: Z_modelo', dp_mod, '   Fe', dp_fe
  print '(a, f14.6, a, f14.6)', 'E_total con apareamiento: Z_modelo', et_mod + dp_mod, '   Fe', et_fe + dp_fe
  print '(a, f14.6)', 'Diferencia con apareamiento (mod - Fe)    =', (et_mod + dp_mod) - (et_fe + dp_fe)
  print '(a, f12.8, a, f12.8)', 'E_B con apareamiento: Z_modelo', eb56 + dp_mod/56.0_dp, '   Fe', eb_fe + dp_fe/56.0_dp

  ! --- E4: fisión simétrica idealizada U-236 -> 2 Pd-118 ---
  q = masa_atomica(92, 236) - 2.0_dp*masa_atomica(46, 118)
  print '(a)', '--- E4 ---'
  print '(a, f16.6)', 'M(92,236)   =', masa_atomica(92, 236)
  print '(a, f16.6)', 'M(46,118)   =', masa_atomica(46, 118)
  print '(a, f16.6)', 'Q [MeV]     =', q
  print '(a, f12.8, a, f12.8)', 'E_B U-236 =', energia_enlace(92, 236), '   E_B Pd-118 =', energia_enlace(46, 118)

contains

  ! ------------------------------------------------------------------
  ! apareamiento: corrección de apareamiento a la energía TOTAL (solo E3)
  !   Entradas: z (protones), a (nucleones)
  !   Salida:   -a_P A^(-1/2) si Z y N pares; +a_P A^(-1/2) si ambos impares;
  !             0 si A es impar [MeV]
  ! ------------------------------------------------------------------
  real(dp) function apareamiento(z, a)
    integer, intent(in) :: z, a
    integer :: n                           ! local: número de neutrones
    n = a - z
    if (mod(a, 2) == 1) then
       apareamiento = 0.0_dp
    else if (mod(z, 2) == 0 .and. mod(n, 2) == 0) then
       apareamiento = -a_p*real(a, dp)**(-0.5_dp)
    else
       apareamiento =  a_p*real(a, dp)**(-0.5_dp)
    end if
  end function apareamiento

end program calculos_E
