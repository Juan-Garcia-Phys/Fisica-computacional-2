! Taller 6 -- Programa principal: fórmula semiempírica de masa
! Usa el módulo modulo_nucleo y escribe tres archivos de resultados:
!   nucleos.dat : Z, A, E_vol, E_sup, E_Coul, E_asim, E_B, M   (8 núclidos)
!   valle.dat   : A, Z de máximo enlace, E_B                   (A = 4, ..., 260)
!   maximo.dat  : A, Z, E_B del núcleo más ligado de valle.dat
program programa_nucleo
  use modulo_nucleo, only: dp, e_volumen, e_superficie, e_coulomb, e_asimetria, &
                           energia_enlace, masa_atomica, z_estable
  implicit none
  integer  :: u_nuc, u_valle, u_max      ! unidades de los tres archivos
  integer  :: a, z                       ! A del ciclo y Z de máximo enlace para ese A
  real(dp) :: eb                         ! E_B del núcleo de máximo enlace para ese A [MeV]
  integer  :: a_max, z_max               ! núcleo más ligado de toda la curva
  real(dp) :: eb_max                     ! su E_B [MeV]

  ! --- 1. Tabla de los ocho núclidos, en el orden del enunciado ---
  open(newunit=u_nuc, file='nucleos.dat', status='replace')
  call escribir_nucleo(u_nuc, 94, 244)   ! Pu-244
  call escribir_nucleo(u_nuc, 32,  70)   ! Ge-70
  call escribir_nucleo(u_nuc, 72, 179)   ! Hf-179
  call escribir_nucleo(u_nuc, 18,  38)   ! Ar-38
  call escribir_nucleo(u_nuc, 26,  56)   ! Fe-56
  call escribir_nucleo(u_nuc, 20,  40)   ! Ca-40
  call escribir_nucleo(u_nuc, 92, 235)   ! U-235
  call escribir_nucleo(u_nuc, 92, 238)   ! U-238
  close(u_nuc)

  ! --- 2. Curva de máximo enlace y búsqueda del núcleo más ligado ---
  ! eb_max empieza en el mayor real posible: el primer E_B lo reemplaza.
  a_max  = 0
  z_max  = 0
  eb_max = huge(1.0_dp)
  open(newunit=u_valle, file='valle.dat', status='replace')
  do a = 4, 260
     call z_estable(a, z, eb)                        ! Z de máximo enlace para este A
     write(u_valle, '(2i6, f18.8)') a, z, eb
     if (eb < eb_max) then                           ! nuevo núcleo más ligado
        a_max  = a
        z_max  = z
        eb_max = eb
     end if
  end do
  close(u_valle)

  ! --- 3. Núcleo más ligado de toda la curva ---
  open(newunit=u_max, file='maximo.dat', status='replace')
  write(u_max, '(2i6, f18.8)') a_max, z_max, eb_max
  close(u_max)

  print '(a)', 'Archivos escritos: nucleos.dat, valle.dat, maximo.dat'
  print '(a, i4, a, i4, a, f12.8, a)', 'Núcleo más ligado: A =', a_max, ', Z =', z_max, &
        ', E_B =', eb_max, ' MeV'

contains

  ! ------------------------------------------------------------------
  ! escribir_nucleo: escribe una fila de nucleos.dat para el núcleo (Z, A)
  !   Entradas: u (unidad del archivo abierto), z (protones), a (nucleones)
  !   Salida:   ninguna; escribe Z, A, los cuatro términos, E_B y M [MeV]
  ! Es subrutina porque realiza una acción (escribir), regla 3.3.
  ! ------------------------------------------------------------------
  subroutine escribir_nucleo(u, z, a)
    integer, intent(in) :: u, z, a
    write(u, '(2i6, 5f18.8, f20.6)') z, a, e_volumen(), e_superficie(a), &
         e_coulomb(z, a), e_asimetria(z, a), energia_enlace(z, a), masa_atomica(z, a)
  end subroutine escribir_nucleo

end program programa_nucleo
