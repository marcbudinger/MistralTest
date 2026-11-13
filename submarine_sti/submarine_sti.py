"""
Traduction Python du modèle Modelica `SubmarineSTI.FullModel`
(chaîne de propulsion sous-marine : moteur DC + arbre + hélice)

Simulation transitoire 0–200 s avec scipy.integrate.solve_ivp (LSODA).

Composants Modelica et leur traduction :
  - TableVoltage           -> v(t) par interpolation linéaire de la table
                              [0,0; 50,700; 100,700; 100,0; 150,0]
                              (palier constant hors table, sémantique Modelica)
  - Resistor R = 8e-3 Ohm, Inductor L = 20e-3 H
  - EMF k = 55.7           : v_emf = k*w ; couple moteur sur l'arbre = k*i
  - Inertia J = 3e3        (moteur)
  - Damper d = 10          (frottement visqueux moteur, lié au bâti)
  - SpringDamper c = 6.5e6 N.m/rad, d = 25000 N.m.s/rad (arbre)
  - BearingFriction [0, 9000] : frottement de Coulomb 9000 N.m
    (signe lissé par tanh pour la robustesse numérique)
  - QuadraticSpeedDependentTorque (MSL 3.2.3, TorqueDirection = true) :
        tau = -tau_nominal*(w/w_nominal)^2*sign(w)
    avec tau_nominal = +1.75e5 N.m à w_nominal = 120 tr/min :
    charge hydrodynamique résistive de l'hélice.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Paramètres (valeurs du fichier .mo)
# ---------------------------------------------------------------------------
R = 8e-3            # Resistor.R              [Ohm]
L = 20e-3           # Inductor.L              [H]
K_EMF = 55.7        # emf.k                   [V.s/rad = N.m/A]
J_MOTOR = 3e3       # inertia_motor.J         [kg.m2]
D_MOTOR = 10.0      # friction_motor.d        [N.m.s/rad]
C_SHAFT = 6.5e6     # shaft_stiffness_loss.c  [N.m/rad]
D_SHAFT = 25000.0   # shaft_stiffness_loss.d  [N.m.s/rad]
J_PROP = 150e3      # Propeler_Inertia.J      [kg.m2]
TAU_BEARING = 9000.0  # Propeler_friction     [N.m]

W_NOMINAL = 120 * 2 * 3.14 / 60   # w_nominal  [rad/s]  (120 tr/min)
TAU_NOMINAL = 1.75e5              # tau_nominal [N.m]  (charge résistive)

T_END = 200.0
W_EPS = 1e-2    # lissage du signe du frottement de Coulomb [rad/s]


def v_source(t):
    """TableVoltage : [0,0; 50,700; 100,700; 100,0; 150,0] (sémantique Modelica)."""
    t = np.asarray(t, dtype=float)
    return np.piecewise(
        t,
        [t <= 0.0, (t > 0.0) & (t < 50.0), (t >= 50.0) & (t < 100.0)],
        [0.0, lambda x: 700.0 * x / 50.0, 700.0, 0.0],
    )


def sgn(w):
    return np.tanh(w / W_EPS)


# ---------------------------------------------------------------------------
# SubmarineSTI.FullModel
# Etats : i [A], w_m [rad/s] (moteur), phi_rel [rad] (torsion d'arbre),
#         w_p [rad/s] (hélice)
# ---------------------------------------------------------------------------
def full_model_rhs(t, y):
    i, w_m, phi_rel, w_p = y
    v = v_source(t)

    di = (v - R * i - K_EMF * w_m) / L

    tau_shaft = C_SHAFT * phi_rel + D_SHAFT * (w_m - w_p)

    dw_m = (K_EMF * i - D_MOTOR * w_m - tau_shaft) / J_MOTOR

    tau_quadratic = -TAU_NOMINAL * (w_p / W_NOMINAL) ** 2 * sgn(w_p)
    tau_friction = -TAU_BEARING * sgn(w_p)
    dw_p = (tau_shaft + tau_friction + tau_quadratic) / J_PROP

    return [di, dw_m, w_m - w_p, dw_p]


def shaft_torque(sol):
    """Couple transmis par l'arbre : c*phi_rel + d*(w_m - w_p)."""
    i, w_m, phi_rel, w_p = sol["y"]
    return C_SHAFT * phi_rel + D_SHAFT * (w_m - w_p)


def simulate_full_model(t_end=T_END, n_points=20001):
    t_eval = np.linspace(0.0, t_end, n_points)
    y0 = [0.0, 0.0, 0.0, 0.0]
    return _simulate_segments(full_model_rhs, y0, t_eval,
                              names=["i", "w_m", "phi_rel", "w_p"])


# ---------------------------------------------------------------------------
# Intégration : découpage aux discontinuités de v(t) (rampe, échelon à t=100)
# ---------------------------------------------------------------------------
def _simulate_segments(rhs, y0, t_eval, names):
    breaks = [0.0, 50.0, 100.0, t_eval[-1]]
    t_all, y_all = [np.array([0.0])], [np.array(y0, dtype=float)]
    y = np.array(y0, dtype=float)
    for a, b in zip(breaks[:-1], breaks[1:]):
        mask = (t_eval > a) & (t_eval <= b)
        if not np.any(mask):
            continue
        te = np.concatenate(([a], t_eval[mask]))
        sol = solve_ivp(rhs, (a, b), y, t_eval=te, method="LSODA",
                        rtol=1e-8, atol=1e-10)
        if not sol.success:
            raise RuntimeError(f"Échec d'intégration sur [{a}, {b}] : {sol.message}")
        t_all.append(sol.t[1:])
        y_all.append(sol.y.T[1:])
        y = sol.y[:, -1]
    t = np.concatenate(t_all)
    Y = np.vstack(y_all)
    return {"t": t, "y": Y.T, "names": names}


def rpm(w):
    return np.asarray(w) * 60.0 / (2 * np.pi)


# ---------------------------------------------------------------------------
# Simulation et figures
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    plt.rcParams.update({"figure.autolayout": True})

    sol = simulate_full_model()
    t = sol["t"]
    i, w_m, phi_rel, w_p = sol["y"]
    tau = shaft_torque(sol)

    # --- Tension et couple transmis par l'arbre ---------------------------
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(t, v_source(t), lw=1.2, color="tab:green")
    ax[0].set_ylabel("Tension [V]")
    ax[0].set_title("FullModel — tension d'alimentation et couple transmis par l'arbre (0–200 s)")
    ax[0].grid(alpha=0.3)
    ax[1].plot(t, tau, lw=1.0, color="tab:orange")
    ax[1].axhline(0.0, color="k", lw=0.5)
    ax[1].set_ylabel("Couple arbre [N.m]")
    ax[1].set_xlabel("Temps [s]")
    ax[1].grid(alpha=0.3)
    fig.savefig("tension_couple_arbre.png", dpi=130)

    # --- Courant et vitesses ----------------------------------------------
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(t, i, lw=1.0)
    ax[0].set_ylabel("Courant i [A]")
    ax[0].set_title("FullModel — courant d'induit et vitesses (0–200 s)")
    ax[0].grid(alpha=0.3)
    ax[1].plot(t, rpm(w_m), lw=1.0, label="moteur")
    ax[1].plot(t, rpm(w_p), "--", lw=1.0, label="hélice")
    ax[1].set_ylabel("Vitesse [tr/min]")
    ax[1].set_xlabel("Temps [s]")
    ax[1].legend()
    ax[1].grid(alpha=0.3)
    fig.savefig("courant_vitesses.png", dpi=130)

    # --- Zoom sur le transitoire de freinage -------------------------------
    fig, ax = plt.subplots(figsize=(10, 5))
    mask = (t >= 95.0) & (t <= 140.0)
    ax.plot(t[mask], tau[mask], lw=1.0, color="tab:orange")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.set_xlabel("Temps [s]")
    ax.set_ylabel("Couple arbre [N.m]")
    ax.set_title("Zoom 95–140 s : transitoire de couple après l'échelon de tension")
    ax.grid(alpha=0.3)
    fig.savefig("couple_arbre_zoom.png", dpi=130)

    idx99 = np.searchsorted(t, 99.0)
    print("Régime établi (t = 99 s) :")
    print(f"  i       = {i[idx99]:9.1f} A")
    print(f"  w_m     = {rpm(w_m[idx99]):8.1f} tr/min")
    print(f"  w_p     = {rpm(w_p[idx99]):8.1f} tr/min")
    print(f"  couple  = {tau[idx99]:9.0f} N.m")
    print(f"  phi_rel = {phi_rel[idx99]:9.4f} rad")
    print(f"Pic de couple après échelon : min = {tau.min():.0f} N.m à t = {t[np.argmin(tau)]:.2f} s, "
          f"max = {tau.max():.0f} N.m à t = {t[np.argmax(tau)]:.2f} s")
    idx_stop = np.argmax((t > 100) & (np.abs(tau) < 100))
    print(f"Retour au repos (|couple| < 100 N.m) : t = {t[idx_stop]:.1f} s")

    plt.show()
