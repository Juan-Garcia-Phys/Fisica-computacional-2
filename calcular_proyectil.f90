%%writefile calcular_proyectil.f90
program calcular_proyectil
    implicit none 

    !declaramos las variables del programa  
    
    real, parameter :: g = 9.8
    real, parameter :: pi = 3.14159265
    real :: v0, angulo_grados, angulo_rad, h_max
    real :: h_prueba, tolerancia
    
    !operaciones de entrada y salida para interactuar con el usuario

    print *, "Ingrese la rapidez inicial en m/s:"
    read *, v0
    print *, "Ingrese el angulo en grados:"
    read *, angulo_grados
    
    !Se cambian los angulos de grados a radianes y se ingresa la formula del mov parabolico
    angulo_rad = angulo_grados * (pi / 180.0) 
    h_max = ((v0**2) * (sin(angulo_rad)**2)) / (2.0 * g)
    
    !interacción de salida
    print *, "La altura maxima es: ", h_max, " m"
    
    !tolerancia que tiene el programa para aprobar el resultado si esta dentro del margen de la milesima
    tolerancia = 0.001 
    
    !calcula nuevamente la altura max pero con el caso de prueba y no con los datos del usuario y se toma el valor absoluto para solo tomar parte positiva
    !9.8m/s angulo de 90
    h_prueba = ((9.8**2) * (sin(90.0 * pi / 180.0)**2)) / (2.0 * g)
    if (abs(h_prueba - 4.9000) <= tolerancia) then
        print *, "Caso 1: PASS"
    else
        print *, "Caso 1: FAIL"
    end if
    !lo mismo pero para el caso dos rapidez de 20.0m/s angulo de 45
    h_prueba = ((20.0**2) * (sin(45.0 * pi / 180.0)**2)) / (2.0 * g)
    if (abs(h_prueba - 10.2041) <= tolerancia) then
        print *, "Caso 2: PASS"
    else
        print *, "Caso 2: FAIL"
    end if
    !lo mismo pero para el caso tres rapidez de 20.0 m/s y un angulo de 30
    h_prueba = ((20.0**2) * (sin(30.0 * pi / 180.0)**2)) / (2.0 * g)
    if (abs(h_prueba - 5.1020) <= tolerancia) then
        print *, "Caso 3: PASS"
    else
        print *, "Caso 3: FAIL"
    end if

end program calcular_proyectil
!gfortran calcular_proyectil.f90 -o calcular_proyectil && ./calcular_proyectil
