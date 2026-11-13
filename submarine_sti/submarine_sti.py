"""
Traduction Python du package Modelica `SubmarineSTI`
(modèles FullModel, ReducedModel1, ReducedModel2 — chaîne de propulsion sous-marine)

Chaque modèle Modelica est réécrit comme un système d'EDO du 1er ordre
intégré avec scipy.integrate.solve_ivp (LSODA).

Composants Modelica et leur traduction :
  - TableVoltage / TimeTable  -> v(t) par interpolation linéaire de la table
                                 [0,0; 50,700; 100,700; 100,0; 150,0]
                                 (extrapolation constante hors table, comme Modelica)
  - Resistor R = 8e-3 Ohm, Inductor L = 20e-3 H
  - EMF k = 55.7  :  v_emf = k*w ; couple moteur sur l'arbre = k*i
  - Inertia J = 3e3 (moteur), J = 150e3 (hélice)
  - Damper d = 10 (frottement visqueux moteur, lié au bâti)
  - SpringDamper c = 6.5e6 N.m/rad, d = 25000 N.m.s/rad (arbre)
  - BearingFriction tau_pos = [0, 9000] : frottement de Coulomb 9000 N.m
    (signe lissé par tanh pour la robustesse numérique, cf. remarque)
  - QuadraticSpeedDependentTorque (MSL 3.2.3) :
        tau = -tau_nominal*(w/w_nominal)^2 * (sign(w) si TorqueDirection sinon 1)
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Paramètres communs (valeurs du fichier .mo)
# ---------------------------------------------------------------------------
R = 8e-3            # Resistor.R              [Ohm]
L = 20e-3           # Inductor.L              [H]
K_EMF = 55.7        # emf.k                   [V.s/rad = N.m/A]
J_MOTOR = 3e3       # inertia_motor.J         [kg.m2]
D_MOTOR = 10.0      # friction_motor.d        [N.m.s/rad]
C_SHAFT = 6.5e6     # shaft_stiffness_loss.c  [N.m/rad]
D_SHAFT = 25000.0   # shaft_stiffness_loss.d  [N.m.s/rad]
J_PROP = 150e3      # Propeler_Inertia.J      [kg.m2]
TAU_BEARING = 9000.0  # Propeler_friction tau_pos = [0, 9000]  [N.m]

W_NOMINAL = 120 * 2 * 3.14 / 60     # w_nominal [rad/s]
TAU_NOMINAL = -1.75e5               # tau_nominal [N.m]

K_RED2 = 700 / (120 * 2 * 3.14 / 60)   # ReducedModel2 emf.k
K_TF = 55.73                           # ReducedModel2 transferFunction

T_END = 150.0

W_EPS = 1e-2   # lissage du signe pour le frottement de Coulomb [rad/s]


def v_source(t):
    """TableVoltage / TimeTable : [0,0; 50,700; 100,700; 100,0; 150,0]
    Interpolation linéaire, palier constant hors table (sémantique Modelica)."""
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

    tau_quadratic = -TAU_NOMINAL * (w_p / W_NOMINAL) ** 2          # TorqueDirection = false
    tau_friction = -TAU_BEARING * sgn(w_p)
    dw_p = (tau_shaft + tau_friction + tau_quadratic) / J_PROP

    return [di, dw_m, w_m - w_p, dw_p]


def simulate_full_model(t_end=T_END, n_points=15001):
    t_eval = np.linspace(0.0, t_end, n_points)
    y0 = [0.0, 0.0, 0.0, 0.0]
    return _simulate_segments(full_model_rhs, y0, t_eval,
                              names=["i", "w_m", "phi_rel", "w_p"])


def shaft_torque(sol):
    """Couple transitoire dans l'arbre : SpringDamper c*phi_rel + d*(w_m - w_p)."""
    i, w_m, phi_rel, w_p = sol["y"]
    return C_SHAFT * phi_rel + D_SHAFT * (w_m - w_p)


# ---------------------------------------------------------------------------
# SubmarineSTI.ReducedModel1
# Etats : i [A], w [rad/s] (inertie hélice rigidement liée à l'EMF)
# (RotationalEMF k = 55.7 ; le heatPort déclaré n'est connecté à rien :
#  transformateur idéal, pertes nulles -> sans effet sur la dynamique)
# ---------------------------------------------------------------------------
def reduced_model1_rhs(t, y):
    i, w = y
    v = v_source(t)

    di = (v - R * i - K_EMF * w) / L

    tau_quadratic = -TAU_NOMINAL * (w / W_NOMINAL) ** 2 * sgn(w)   # TorqueDirection = true
    tau_friction = -TAU_BEARING * sgn(w)
    dw = (K_EMF * i - D_MOTOR * w + tau_friction + tau_quadratic) / J_PROP

    return [di, dw]


def simulate_reduced_model1(t_end=T_END, n_points=15001):
    t_eval = np.linspace(0.0, t_end, n_points)
    y0 = [0.0, 0.0]
    return _simulate_segments(reduced_model1_rhs, y0, t_eval,
                              names=["i", "w"])


# ---------------------------------------------------------------------------
# SubmarineSTI.ReducedModel2
# Partie électrique seule (i, w) + bloc TransferFunction :
#   num b = {0, J/k, 0}, den a = {L*J/k^2, R*J/k^2, 1},  k = 55.73
#   G(s) = (J/k)*s / ( (L*J/k^2)*s^2 + (R*J/k^2)*s + 1 )   -> courant i(t)
# ---------------------------------------------------------------------------
A2_TF = L * J_PROP / K_TF ** 2
A1_TF = J_PROP * R / K_TF ** 2
B1_TF = J_PROP / K_TF


def reduced_model2_rhs(t, y):
    i, w, x1, x2 = y
    v = v_source(t)

    di = (v - R * i - K_RED2 * w) / L
    dw = (K_RED2 * i) / J_PROP

    dx1 = x2
    dx2 = (v - x1 - A1_TF * x2) / A2_TF
    y_tf = B1_TF * x2

    return [di, dw, dx1, dx2], y_tf


def simulate_reduced_model2(t_end=T_END, n_points=15001):
    t_eval = np.linspace(0.0, t_end, n_points)

    def rhs(t, y):
        return reduced_model2_rhs(t, y)[0]

    y0 = [0.0, 0.0, 0.0, 0.0]
    sol = _simulate_segments(rhs, y0, t_eval,
                             names=["i", "w", "x1", "x2"])
    sol["y_tf"] = np.array([reduced_model2_rhs(t, y)[1]
                            for t, y in zip(sol["t"], sol["y"].T)])
    return sol


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
    return w * 60.0 / (2 * np.pi)


# ---------------------------------------------------------------------------
# Exécution des trois simulations transitoires et tracés
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    plt.rcParams.update({"figure.autolayout": True})

    # --- FullModel (700 s) -------------------------------------------------
    s1 = simulate_full_model(t_end=700.0, n_points=70001)
    t = s1["t"]
    i, w_m, phi_rel, w_p = s1["y"]

    tau_shaft = shaft_torque(s1)

    fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    ax[0].plot(t, i, lw=1.2)
    ax[0].set_ylabel("Courant i [A]")
    ax[0].set_title("FullModel — transitoire (0–150 s)")
    ax[1].plot(t, rpm(w_m), lw=1.2, label="moteur")
    ax[1].plot(t, rpm(w_p), lw=1.2, label="hélice")
    ax[1].set_ylabel("Vitesse [tr/min]")
    ax[1].legend()
    ax[2].plot(t, v_source(t), lw=1.2, label="v(t) [V]")
    ax_t = ax[2].twinx()
    ax_t.plot(t, tau_shaft, "tab:orange", lw=1.2, label="couple arbre [N.m]")
    ax_t.set_ylabel("Couple arbre [N.m]")
    ax[2].set_ylabel("Tension [V]")
    ax[2].set_xlabel("Temps [s]")
    fig.savefig("full_model.png", dpi=130)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, tau_shaft, lw=1.0, label="couple total arbre")
    ax.plot(t, C_SHAFT * phi_rel, "--", lw=1.0, label="partie ressort c*phi_rel")
    ax.plot(t, D_SHAFT * (w_m - w_p), ":", lw=1.0, label="partie amortissement d*(w_m - w_p)")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.axvline(100.0, color="gray", lw=0.8, alpha=0.6)
    ax.set_xlabel("Temps [s]")
    ax.set_ylabel("Couple [N.m]")
    ax.set_title("FullModel — couple transitoire dans l'arbre (0–700 s)")
    ax.legend()
    fig.savefig("shaft_torque_700s.png", dpi=130)

    for ts in (50.0, 99.0, 100.0, 105.0, 150.0, 300.0, 700.0):
        idx = np.searchsorted(t, min(ts, t[-1]))
        print(f"  t = {ts:6.1f} s : couple arbre = {tau_shaft[idx]:12.1f} N.m")

    print("FullModel — régime établi (t=99 s) :")
    print(f"  i      = {i[np.searchsorted(t, 99)]:10.1f} A")
    print(f"  w_m    = {rpm(w_m[np.searchsorted(t, 99)]):10.1f} tr/min")
    print(f"  w_p    = {rpm(w_p[np.searchsorted(t, 99)]):10.1f} tr/min")
    print(f"  phi_rel= {phi_rel[np.searchsorted(t, 99)]:10.4f} rad")

    s1b = simulate_full_model()
    t99 = np.searchsorted(s1b["t"], 99)
    print(f"  i(99s) = {s1b['y'][0][t99]:10.1f} A")

    # --- ReducedModel1 -----------------------------------------------------
    s2 = simulate_reduced_model1()
    t2 = s2["t"]
    i2, w2 = s2["y"]

    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax[0].plot(t2, i2, lw=1.2)
    ax[0].set_ylabel("Courant i [A]")
    ax[0].set_title("ReducedModel1 — transitoire (0–150 s)")
    ax[1].plot(t2, rpm(w2), lw=1.2)
    ax[1].set_ylabel("Vitesse [tr/min]")
    ax[1].set_xlabel("Temps [s]")
    fig.savefig("reduced_model1.png", dpi=130)

    print("ReducedModel1 — régime établi (t=99 s) :")
    print(f"  i = {i2[np.searchsorted(t2, 99)]:10.1f} A")
    print(f"  w = {rpm(w2[np.searchsorted(t2, 99)]):10.1f} tr/min")

    # --- ReducedModel2 -----------------------------------------------------
    s3 = simulate_reduced_model2()
    t3 = s3["t"]
    i3, w3 = s3["y"][:2]
    y_tf = s3["y_tf"]

    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax[0].plot(t3, i3, lw=1.2, label="circuit électrique")
    ax[0].plot(t3, y_tf, "--", lw=1.2, label="TransferFunction")
    ax[0].set_ylabel("Courant i [A]")
    ax[0].set_title("ReducedModel2 — transitoire (0–150 s)")
    ax[0].legend()
    ax[1].plot(t3, rpm(w3), lw=1.2)
    ax[1].set_ylabel("Vitesse [tr/min]")
    ax[1].set_xlabel("Temps [s]")
    fig.savefig("reduced_model2.png", dpi=130)

    print("ReducedModel2 — régime établi (t=99 s) :")
    print(f"  i (circuit)        = {i3[np.searchsorted(t3, 99)]:10.1f} A")
    print(f"  i (TransferFunction)= {y_tf[np.searchsorted(t3, 99)]:10.1f} A")
    print(f"  w                  = {rpm(w3[np.searchsorted(t3, 99)]):10.1f} tr/min")

    plt.show()
