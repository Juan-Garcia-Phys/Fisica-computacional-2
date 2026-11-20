# ============================================================
# Programa:  maximo_espectro.py   (Taller 8, Parte D)
# Propósito: hallar la longitud de onda del máximo del espectro de
#            un cuerpo negro a T = 5772 K (modelo del Sol) en mallas
#            de paso 100, 50, 25 y 12.5 nm, de tres maneras, y medir
#            el error frente a la ley de Wien, lambda_max = b/T.
#   (a) mayor valor del arreglo
#   (b) cero del spline de la derivada central (derivada programada)
#   (c) cero de la derivada del spline de B_lambda
# Entrada:   constante_wien.dat (x5, x3, b) de la Parte B
# Salida:    maximo_espectro.dat, espectro_solar.png
# Restricción del taller: no se usa np.gradient ni scipy.optimize;
#            CubicSpline (.roots y .derivative) solo en esta parte.
# ============================================================
import numpy as np
import matplotlib
matplotlib.use("Agg")                    # guarda la figura sin ventana
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Constantes exactas del SI
H, C, KB = 6.62607015e-34, 299792458.0, 1.380649e-23
T_SOL = 5772.0                           # K
PASOS = [100.0, 50.0, 25.0, 12.5]        # nm, en el orden que pide el enunciado
AZUL, DORADO, ROJO, GRIS = "#162A54", "#C9972C", "#B02A20", "#7F7F7F"


def planck_lambda(lam_nm, T):
    """Radiancia espectral B_lambda [W m^-2 sr^-1 m^-1] con lambda en nm."""
    lam = lam_nm*1e-9                                  # nm -> m
    x = H*C/(lam*KB*T)                                 # energía del fotón / kB T
    return 2*H*C**2/lam**5/np.expm1(x)                 # expm1 = e^x - 1 sin cancelación


def derivada_central(y, h):
    """Derivada central de orden 2 en los puntos interiores (sección 5.2)."""
    return (y[2:] - y[:-2])/(2*h)


def elegir_maximo(spline_derivada, referencia):
    """Ceros de la derivada con segunda derivada negativa (máximos).

    En una malla gruesa el spline oscila y produce ceros espurios; si
    queda más de un máximo se conserva el más cercano a 'referencia',
    el resultado del método (a).
    """
    r = spline_derivada.roots(extrapolate=False)
    r = r[np.isfinite(r)]
    r = r[spline_derivada(r, 1) < 0]                   # f'' < 0: máximo
    if len(r) == 0:
        raise RuntimeError("no se encontró un máximo en la malla")
    return r[np.argmin(np.abs(r - referencia))]


def maximos(paso):
    """Las tres estimaciones de lambda_max [nm] para una malla de paso dado."""
    lam = np.arange(100.0, 3000.0 + paso/2, paso)      # + paso/2 incluye 3000 nm
    B = planck_lambda(lam, T_SOL)
    lam_a = lam[np.argmax(B)]                                          # (a)
    s_d = CubicSpline(lam[1:-1], derivada_central(B, paso))            # (b)
    lam_b = elegir_maximo(s_d, lam_a)
    dB = CubicSpline(lam, B).derivative()                              # (c)
    lam_c = elegir_maximo(dB, lam_a)
    return lam, B, (lam_a, lam_b, lam_c)


# --- Constantes de la Parte B (Fortran)
x5, x3, b_wien = np.loadtxt("constante_wien.dat")
lam_ref = b_wien/T_SOL*1e9                             # ley de Wien, nm
lam_nu = H*C/(KB*T_SOL*x3)*1e9                         # máximo de B_nu, expresado en lambda
print(f"lambda_max = b/T = {lam_ref:.4f} nm;  máximo en frecuencia -> lambda = {lam_nu:.1f} nm\n")

# --- Errores absolutos en las cuatro mallas
filas = []
print(" paso[nm]  lambda(a)    error(a)   lambda(b)    error(b)   lambda(c)    error(c)")
for paso in PASOS:
    lam, B, estimaciones = maximos(paso)
    errores = [abs(e - lam_ref) for e in estimaciones]
    filas.append([paso, *errores])
    print(f"{paso:8.1f}" + "".join(f"{e:11.4f}{err:12.3e}" for e, err in zip(estimaciones, errores)))
    if paso == 25.0:
        lam25, B25, est25 = lam, B, estimaciones       # malla de la figura

filas = np.array(filas)
np.savetxt("maximo_espectro.dat", filas, fmt="%24.16e")
print("\nrazón entre errores consecutivos (b) y (c):")
print(np.round(filas[:-1, 2:]/filas[1:, 2:], 2))

# --- Figura: espectro con la malla de 25 nm
fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=130)
lam_fino = np.linspace(100, 3000, 2000)
ax.plot(lam_fino, planck_lambda(lam_fino, T_SOL)*1e-9, color=GRIS, lw=2.5, alpha=0.5,
        label="Planck exacta")
ax.plot(lam25, B25*1e-9, "o", ms=3, color=AZUL, label=r"malla $\Delta\lambda = 25$ nm")
ax.axvline(lam_ref, color=DORADO, ls="--",
           label=rf"$\lambda_{{\max}} = b/T = {lam_ref:.1f}$ nm")
ax.plot(est25[2], planck_lambda(est25[2], T_SOL)*1e-9, "*", ms=12, color=DORADO,
        label=rf"máximo con spline (c): {est25[2]:.1f} nm")
ax.axvline(lam_nu, color=ROJO, ls=":",
           label=rf"máximo en frecuencia: $hc/(k_BTx_3) = {lam_nu:.0f}$ nm")
ax.set_xlabel(r"$\lambda$ [nm]")
ax.set_ylabel(r"$B_\lambda$ [W m$^{-2}$ sr$^{-1}$ nm$^{-1}$]")
ax.set_title(r"Espectro de cuerpo negro a $T = 5772$ K (Sol)")
ax.set_xlim(100, 3000); ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("espectro_solar.png")
print("\nfigura guardada en espectro_solar.png")
