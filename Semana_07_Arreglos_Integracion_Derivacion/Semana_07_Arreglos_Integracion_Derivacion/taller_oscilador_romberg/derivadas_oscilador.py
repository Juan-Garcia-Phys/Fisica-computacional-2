# ============================================================
# Programa:  derivadas_oscilador.py   (Taller 7, Parte C)
# Propósito: velocidad y aceleración del oscilador amortiguado a
#            partir de las posiciones, con diferencias finitas
#            sobre el arreglo completo (rebanadas, sin ciclos sobre
#            los puntos), comparadas con las derivadas exactas.
# Salida:    derivadas_oscilador.dat  (dt = 0.01: t, x, v, a centrales)
#            errores_derivadas.dat    (dt y cinco errores máximos)
#            derivadas_oscilador.png  (v y a numéricas y exactas)
# ============================================================
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAMMA = 0.5                               # amortiguamiento [1/s]
OMEGA = 2*np.pi                           # frecuencia angular [rad/s]
T_FINAL = 3.0                             # [s]


def posicion(t):
    """x(t) = exp(-gamma t) cos(omega t) [m]."""
    return np.exp(-GAMMA*t)*np.cos(OMEGA*t)


def velocidad_exacta(t):
    """v(t) exacta [m/s], solo para comparar."""
    return -np.exp(-GAMMA*t)*(GAMMA*np.cos(OMEGA*t) + OMEGA*np.sin(OMEGA*t))


def aceleracion_exacta(t):
    """a(t) exacta [m/s^2], solo para comparar."""
    return np.exp(-GAMMA*t)*((GAMMA**2 - OMEGA**2)*np.cos(OMEGA*t)
                             + 2*GAMMA*OMEGA*np.sin(OMEGA*t))


def malla(dt):
    """Tiempos en [0, T_FINAL] con paso dt y posiciones x(t)."""
    n = int(round(T_FINAL/dt)) + 1
    t = dt*np.arange(n)
    return t, posicion(t)


def derivadas(x, dt):
    """Diferencias finitas sobre todo el arreglo. Devuelve cada resultado
    con la rebanada de t donde está definido:
      adelantada:   índices 0..n-2   necesita x[i+1]
      atrasada:     índices 1..n-1   necesita x[i-1]
      central (v,a): índices 1..n-2   necesita x[i-1] y x[i+1]
      cinco puntos: índices 2..n-3   necesita x[i-2] ... x[i+2]"""
    adel = (x[1:] - x[:-1])/dt
    atras = (x[1:] - x[:-1])/dt          # mismos cocientes, asignados al punto derecho
    cent = (x[2:] - x[:-2])/(2*dt)
    cinco = (-x[4:] + 8*x[3:-1] - 8*x[1:-3] + x[:-4])/(12*dt)
    a_cent = (x[2:] - 2*x[1:-1] + x[:-2])/dt**2
    return {
        "adelantada": (adel, slice(0, -1)),
        "atrasada": (atras, slice(1, None)),
        "central": (cent, slice(1, -1)),
        "cinco_puntos": (cinco, slice(2, -2)),
        "a_central": (a_cent, slice(1, -1)),
    }


def errores_maximos(dt):
    """Error máximo de cada fórmula frente a la derivada exacta."""
    t, x = malla(dt)
    d = derivadas(x, dt)
    errores = []
    for nombre in ["adelantada", "atrasada", "central", "cinco_puntos"]:
        valores, idx = d[nombre]
        errores.append(np.max(np.abs(valores - velocidad_exacta(t[idx]))))
    valores, idx = d["a_central"]
    errores.append(np.max(np.abs(valores - aceleracion_exacta(t[idx]))))
    return errores


# --- 1. dt = 0.01: archivo de datos y figura con las fórmulas centrales
dt = 0.01
t, x = malla(dt)
d = derivadas(x, dt)
v_c, idx = d["central"]
a_c, _ = d["a_central"]
t_int, x_int = t[idx], x[idx]                       # 299 puntos interiores
np.savetxt("derivadas_oscilador.dat", np.column_stack([t_int, x_int, v_c, a_c]),
           fmt="%.15e", header="t [s]   x [m]   v central [m/s]   a central [m/s^2]  (dt = 0.01 s)")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
t_fino = np.linspace(0, T_FINAL, 2000)
ax1.plot(t_fino, velocidad_exacta(t_fino), color="#162A54", lw=1.5, label="exacta")
ax1.plot(t_int[::5], v_c[::5], "o", ms=3, color="#C9972C", label="central, $\\Delta t$ = 0.01 s")
ax1.set_ylabel("$v$ [m/s]"); ax1.legend(loc="upper right"); ax1.grid(alpha=0.3)
ax2.plot(t_fino, aceleracion_exacta(t_fino), color="#162A54", lw=1.5, label="exacta")
ax2.plot(t_int[::5], a_c[::5], "o", ms=3, color="#B02A20", label="central, $\\Delta t$ = 0.01 s")
ax2.set_ylabel("$a$ [m/s$^2$]"); ax2.set_xlabel("$t$ [s]"); ax2.legend(loc="upper right"); ax2.grid(alpha=0.3)
fig.suptitle("Oscilador amortiguado: derivadas numéricas frente a exactas")
fig.tight_layout()
fig.savefig("derivadas_oscilador.png", dpi=150)

# --- 2. Errores máximos con dos pasos y razones (orden de cada fórmula)
pasos = [0.01, 0.005]
tabla = np.array([[p] + errores_maximos(p) for p in pasos])
np.savetxt("errores_derivadas.dat", tabla, fmt="%.6e",
           header="dt  err_adelantada  err_atrasada  err_central  err_cinco_puntos  err_a_central")

nombres = ["v adelantada", "v atrasada", "v central", "v cinco puntos", "a central"]
esperadas = [2, 2, 4, 16, 4]
print(f"{'fórmula':>15} {'error dt=0.01':>14} {'error dt=0.005':>15} {'razón':>7} {'esperada':>9}")
for i, nombre in enumerate(nombres):
    e1, e2 = tabla[0, i + 1], tabla[1, i + 1]
    print(f"{nombre:>15} {e1:14.4e} {e2:15.4e} {e1/e2:7.2f} {esperadas[i]:9d}")
