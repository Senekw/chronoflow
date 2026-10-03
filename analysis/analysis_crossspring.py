"""Cross-spring flexure pivot under gravity load (design alternative that was REJECTED).
Euler-Bernoulli beam-column FE (elastic + geometric stiffness) for each leaf; the balance
rotates rigidly about the crossing point C located a fraction lam of the leaf length from
the fixed clamp.  k_theta(P) = 2 U(theta=1).  P > 0 = axial compression."""
import numpy as np

def beam_column_stiffness(EI, L, P, n_el=60):
    le = L / n_el
    ke = (EI / le**3) * np.array([[12, 6*le, -12, 6*le], [6*le, 4*le**2, -6*le, 2*le**2],
                                  [-12, -6*le, 12, -6*le], [6*le, 2*le**2, -6*le, 4*le**2]])
    kg = (1 / (30 * le)) * np.array([[36, 3*le, -36, 3*le], [3*le, 4*le**2, -3*le, -le**2],
                                     [-36, -3*le, 36, -3*le], [3*le, -le**2, -3*le, 4*le**2]])
    n = 2 * (n_el + 1); K = np.zeros((n, n))
    for e in range(n_el):
        idx = [2*e, 2*e+1, 2*e+2, 2*e+3]; K[np.ix_(idx, idx)] += ke - P * kg
    return K

def leaf_energy(EI, L, P, lam, theta=1.0, n_el=60):
    K = beam_column_stiffness(EI, L, P, n_el); n = K.shape[0]
    d = np.zeros(n); d[-2] = theta * (L - lam * L); d[-1] = theta
    fixed = [0, 1, n-2, n-1]; free = [i for i in range(n) if i not in fixed]
    d[free] = np.linalg.solve(K[np.ix_(free, free)], -K[np.ix_(free, fixed)] @ d[fixed])
    return 0.5 * d @ K @ d

def pivot_stiffness(EI, L, lam, loads=(0.0, 0.0)):
    return sum(2 * leaf_energy(EI, L, P, lam) for P in loads)

def gravity_sensitivity(E, w, t, L, m_bal, lam=0.5):
    """Worst-case fractional stiffness change over all in-plane gravity directions."""
    EI = E * w * t**3 / 12; k0 = pivot_stiffness(EI, L, lam); worst = 0.0
    e1 = np.array([np.cos(np.pi/4), np.sin(np.pi/4)]); e2 = np.array([np.cos(-np.pi/4), np.sin(-np.pi/4)])
    for ang in np.deg2rad(np.arange(0, 360, 15)):
        g = m_bal * 9.81 * np.array([np.cos(ang), np.sin(ang)])
        kl = pivot_stiffness(EI, L, lam, (-(g @ e1), -(g @ e2)))
        worst = max(worst, abs(kl / k0 - 1))
    return worst, k0, EI
