"""Taller 6 -- Módulo de análisis de la curva de máximo enlace.

Trabaja solo con los archivos que escribe el programa de Fortran (no recalcula
la fórmula semiempírica). Convención: E_B < 0 para un núcleo ligado; en la
gráfica se usa B/A = -E_B, positiva, como en las tablas experimentales.
"""
import numpy as np
import matplotlib.pyplot as plt

AZUL, DORADO, ROJO = "#162A54", "#C9972C", "#B02A20"     # colores del curso


def leer_valle(archivo):
    """Lee valle.dat, escrito por programa_nucleo.f90.

    Entrada: archivo -- ruta de valle.dat
    Salida:  tupla (A, Z, E_B): A y Z como enteros, E_B [MeV] como reales;
             cada una es una columna del archivo (un arreglo de NumPy)
    """
    datos = np.loadtxt(archivo)
    a = datos[:, 0].astype(int)
    z = datos[:, 1].astype(int)
    e_b = datos[:, 2]
    return a, z, e_b


def nucleo_mas_ligado(a, z, e_b):
    """Encuentra el núcleo con mayor enlace por nucleón (E_B más negativo).

    Recorre los datos con un ciclo y un condicional, como en Fortran
    (sin min, np.min ni np.argmin).
    Entradas: columnas A, Z y E_B [MeV] de leer_valle
    Salida:   tupla (A, Z, E_B) del núcleo más ligado
    """
    # El primer núcleo es el mínimo provisional; los demás compiten con él
    a_max, z_max, eb_max = int(a[0]), int(z[0]), float(e_b[0])
    for a_i, z_i, eb_i in zip(a[1:], z[1:], e_b[1:]):
        if eb_i < eb_max:                       # más negativo = más ligado
            a_max, z_max, eb_max = int(a_i), int(z_i), float(eb_i)
    return a_max, z_max, eb_max


def escribir_maximo(archivo, a_max, z_max, eb_max):
    """Escribe el núcleo más ligado con el mismo formato de maximo.dat.

    Entradas: archivo -- ruta de salida (maximo_python.dat);
              a_max, z_max, eb_max -- A, Z y E_B [MeV] del núcleo
    Salida:   ninguna; escribe una fila "A Z E_B" (formato Fortran 2i6, f18.8)
    """
    with open(archivo, "w") as f:
        f.write(f"{a_max:6d}{z_max:6d}{eb_max:18.8f}\n")


def graficar_curva(a, e_b, experimentales, archivo):
    """Grafica B/A = -E_B contra A con los valores experimentales y la guarda.

    Entradas: a, e_b -- columnas A y E_B [MeV] de la curva de máximo enlace;
              experimentales -- lista de tuplas (símbolo, Z, A, B/A [MeV]);
              archivo -- ruta de la imagen (curva_enlace.png)
    Salida:   la figura de matplotlib (también queda guardada en archivo)
    """
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(a, -e_b, "-", color=AZUL, lw=2,
            label="Modelo de 4 términos (máximo enlace para cada $A$)")
    # Puntos experimentales: cada uno con su propio Z (comparación orientativa)
    a_exp = [fila[2] for fila in experimentales]
    ba_exp = [fila[3] for fila in experimentales]
    ax.plot(a_exp, ba_exp, "o", color=DORADO, mec=ROJO, ms=7,
            label="Experimental (AME2020)")
    # Desplazamientos de las etiquetas que quedarían encimadas (Fe/Ni y U-235/U-238)
    desplazamientos = {"Fe-56": (-30, 6), "Ni-62": (4, 6), "U-235": (-14, 8), "U-238": (4, -12)}
    for simbolo, _, a_i, ba_i in experimentales:
        ax.annotate(simbolo, (a_i, ba_i), textcoords="offset points",
                    xytext=desplazamientos.get(simbolo, (5, -12)), fontsize=8)
    ax.set_xlabel("Número de nucleones $A$")
    ax.set_ylabel("Energía de enlace por nucleón $B/A = -E_B$ [MeV]")
    ax.set_title("Curva de máximo enlace: fórmula semiempírica frente a datos")
    ax.set_ylim(6.0, 9.1)                        # zona de interés; el modelo cae fuera para A muy pequeño
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(archivo, dpi=150)
    return fig
