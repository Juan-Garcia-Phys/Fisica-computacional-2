! Ajuste lineal por minimos cuadrados de T^2 vs L para estimar g (pendulo simple)
program ajuste_pendulo
  use iso_fortran_env, only: real64, iostat_end   ! Doble precisión y código de fin de archivo
  implicit none                                   ! Obliga a declarar todas las variables

  integer :: unidad, estado, unidad_salida        ! Unidades de archivo y estado de lectura
  integer :: id, n_osc, n                         ! Campos enteros de cada fila y contador N
  real(real64) :: long_cm, angulo, tiempo_s       ! Campos reales leídos de cada fila
  real(real64) :: L, T, x, y                      ! Longitud (m), período (s), x = L, y = T^2
  real(real64) :: sx, sy, sxx, sxy                ! Acumuladores de mínimos cuadrados
  real(real64) :: a, b, g, pi                     ! Pendiente, intercepto, gravedad y pi
  real(real64) :: den, y_media, y_ajustada, sse, sst, r2  ! Variables para la calidad del ajuste

  pi = 4.0_real64 * atan(1.0_real64)              ! pi con precisión doble: atan(1) = pi/4

  ! Todo acumulador debe empezar en cero antes del ciclo
  n = 0
  sx = 0.0_real64
  sy = 0.0_real64
  sxx = 0.0_real64
  sxy = 0.0_real64

  ! Abrir el archivo limpio generado
  open(newunit=unidad, file='pendulo_limpio.dat', status='old', &
       action='read', iostat=estado)
  if (estado /= 0) error stop 'Error: No fue posible abrir pendulo_limpio.dat'

  ! Primera lectura: calcular sumatorias
  ! El número de filas es desconocido: se lee hasta encontrar el final del archivo
  do
    read(unidad, *, iostat=estado) id, long_cm, angulo, n_osc, tiempo_s
    if (estado == iostat_end) exit                ! Fin normal del archivo
    if (estado /= 0) error stop 'Error de lectura' ! Cualquier otro valor indica una fila mal formada

    L = long_cm / 100.0_real64                    ! Convierte cm a m: g debe salir en m/s^2
    T = tiempo_s / real(n_osc, real64)            ! Período de UNA oscilación (N puede ser 5, 10, 15 o 20)
    x = L                                         ! Variable independiente de la recta
    y = T**2                                      ! Linealización: T^2 = (4 pi^2 / g) L

    n = n + 1                                     ! Cuenta la fila aceptada
    sx = sx + x
    sy = sy + y
    sxx = sxx + x**2
    sxy = sxy + x * y
  end do

  ! Con menos de dos puntos no se puede trazar una recta
  if (n < 2) error stop 'Se requieren al menos dos datos'

  ! Si todas las longitudes son iguales, el denominador es cero y la pendiente no existe
  den = real(n, real64)*sxx - sx*sx
  if (abs(den) <= tiny(den)) error stop 'Denominador nulo'

  a = (real(n, real64)*sxy - sx*sy) / den         ! Pendiente (s^2/m)
  b = (sy - a*sx) / real(n, real64)               ! Intercepto (s^2); el modelo ideal predice b = 0
  g = (4.0_real64 * pi**2) / a                    ! Despeja g de a = 4 pi^2 / g

  ! Segunda lectura: calcular R^2
  ! R^2 necesita la media de y y la recta ya calculadas, por eso exige recorrer el archivo otra vez
  rewind(unidad)                                  ! Vuelve al inicio del archivo sin cerrarlo
  sse = 0.0_real64
  sst = 0.0_real64
  y_media = sy / real(n, real64)

  do
    read(unidad, *, iostat=estado) id, long_cm, angulo, n_osc, tiempo_s
    if (estado == iostat_end) exit
    if (estado /= 0) error stop 'Error de lectura en la segunda pasada'  ! Evita un R^2 erroneo sin aviso

    L = long_cm / 100.0_real64                    ! Mismas conversiones que en la primera lectura
    T = tiempo_s / real(n_osc, real64)
    x = L
    y = T**2

    y_ajustada = a * x + b                        ! Predicción de la recta
    sse = sse + (y - y_ajustada)**2               ! Suma de residuos al cuadrado
    sst = sst + (y - y_media)**2                  ! Variación total alrededor de la media
  end do
  close(unidad)                                   ! Libera el archivo de datos

  ! R^2 solo está definido si los valores de y no son todos iguales
  if (sst > 0.0_real64) then
    r2 = 1.0_real64 - (sse / sst)
  else
    r2 = 0.0_real64
  end if


  print '(A)',          '--- RESULTADOS DEL AJUSTE EN FORTRAN ---'
  print '(A,I0)',       'N                     = ', n
  print '(A,F10.6,A)',  'Pendiente a           = ', a, ' s^2/m'
  print '(A,F10.6,A)',  'Intercepto b          = ', b, ' s^2'
  print '(A,F10.6)',    'R^2                   = ', r2
  print '(A,F10.4,A)',  'Gravedad calculada g  = ', g, ' m/s^2'
  print '(A)',          '----------------------------------------'

  ! Guardar los resultados: cinco valores, en este orden y sin encabezado
  open(newunit=unidad_salida, file='resultados_ajuste.dat', &
       status='replace', action='write')
  write(unidad_salida, *) n, a, b, r2, g
  close(unidad_salida)

end program ajuste_pendulo
