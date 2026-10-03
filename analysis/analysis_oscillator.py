"""Taut-band balance oscillator analysis for ChronoFlow.

The balance wheel is clamped at the midpoint of a thin spring-steel ribbon (feeler-gauge
stock) that is clamped to the frame at both ends and pre-tensioned by a steel coil spring.
Each half-band is a torsion spring:

    k_half(T) = G*J/L_h + T*w^2/(12*L_h),   J = w*t^3/3 * (1 - 0.63 t/w)     (Timoshenko)

Total k = k_upper + k_lower, which depends on the SUM of the two half-band tensions.
The tensioner spring fixes the tension of its own half at T0; with the band vertical the
other (fixed-clamp) half carries T0 +/- m g.  The gravity effect is therefore FIRST order,
Delta k = +/- m g w^2 / (12 L_h), but small because the tension term is only ~9 % of k and
m g << T0 (computed: +/-0.16 % in frequency relative to the band-horizontal case, i.e.
0.32 % peak-to-peak between fixed-clamp-up and fixed-clamp-down).  There is no pivot
friction and no lubricant, and the restoring torque is linear -> intrinsically isochronous.
"""
import numpy as np

G_STEEL = 80e9          # Pa shear modulus, spring steel
E_STEEL = 200e9
ALPHA_G_STEEL = -2.6e-4  # 1/K, temperature coefficient of G
ALPHA_PETG = 6.0e-5      # 1/K linear expansion (balance ring)
TAU_FATIGUE_STEEL = 350e6  # Pa, conservative fully-reversed torsional fatigue limit for hardened shim


def band_J(w, t):
    return w * t**3 / 3 * (1 - 0.63 * t / w)


def k_half(w, t, L_h, T, G=G_STEEL):
    return G * band_J(w, t) / L_h + T * w**2 / (12 * L_h)


def k_total(w, t, L_h, T_upper, T_lower, G=G_STEEL):
    return k_half(w, t, L_h, T_upper, G) + k_half(w, t, L_h, T_lower, G)


def ring_inertia(r_out, r_in, thick, rho):
    m = rho * np.pi * (r_out**2 - r_in**2) * thick
    I = 0.5 * m * (r_out**2 + r_in**2)
    return m, I


def balance_inertia(r_out, r_in, thick, rho, n_trim, m_trim, r_trim, hub_m=0.0, hub_r=0.0):
    m_ring, I_ring = ring_inertia(r_out, r_in, thick, rho)
    I_trim = n_trim * m_trim * r_trim**2
    I_hub = 0.5 * hub_m * hub_r**2
    return m_ring + n_trim * m_trim + hub_m, I_ring + I_trim + I_hub


def frequency(k, I):
    return np.sqrt(k / I) / (2 * np.pi)


def max_shear_stress(t, theta, L_h, G=G_STEEL):
    """Peak shear stress in a thin ribbon twisted by theta over length L_h."""
    return G * t * theta / L_h


def orientation_study(w, t, L_h, T0, m_bal, I_bal):
    """Frequency in the three principal orientations.
    Tensioner spring on the LOWER end in 'upright'. Band axis vertical: the half between
    balance and the fixed clamp carries T0 + m g (upright) or T0 - m g (inverted);
    the tensioner-side half carries T0 always.  Band horizontal: both halves T0 (gravity
    only sags the band)."""
    g = 9.81
    cases = {
        "band horizontal (device flat)": (T0, T0),
        "band vertical, fixed clamp up": (T0 + m_bal * g, T0),
        "band vertical, fixed clamp down": (T0 - m_bal * g, T0),
    }
    out = {}
    for name, (Tu, Tl) in cases.items():
        k = k_total(w, t, L_h, Tu, Tl)
        out[name] = frequency(k, I_bal)
    return out


def band_sag(m_bal, L_h, T0):
    """Lateral sag of the balance when the band is horizontal (point load at mid-span)."""
    return m_bal * 9.81 * L_h / (2 * 2 * T0)  # two halves share the load: F*L_h/(4 T) ... conservative


def temperature_coefficient(w, t, L_h, T0, alpha_tensioner=0.0):
    """d(ln f)/dT [1/K]: shear modulus drift + ring expansion + tensioner drift (steel spring ~0)."""
    kG = 2 * G_STEEL * band_J(w, t) / L_h
    kT = 2 * T0 * w**2 / (12 * L_h)
    k = kG + kT
    dk_k = (kG * ALPHA_G_STEEL + kT * alpha_tensioner) / k
    dI_I = 2 * ALPHA_PETG
    return 0.5 * (dk_k - dI_I)


if __name__ == "__main__":
    import design_params as dp
    k = k_total(dp.BAND_W, dp.BAND_T, dp.BAND_LH, dp.BAND_T0, dp.BAND_T0)
    print(f"k_total = {k:.4e} N m/rad")
