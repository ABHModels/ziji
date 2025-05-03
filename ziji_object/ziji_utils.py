import numpy as np
import math
from scipy.integrate import quad
from math import pi, sqrt
import warnings
from scipy.integrate import cumtrapz


def BlackBodyIntensity(energy_keV, kbTeff_keV, f_color_factor, Y_angle_func):
    # Assuming constant values are defined somewhere
    h_plank_cgs = 6.62607015e-27  # Planck constant in cgs units
    c_light_cgs = 2.99792458e10  # Speed of light in cgs units
    keV_erg = 1.60218e-9  # 1 keV in ergs

    A1 = 2.0 / f_color_factor**4 * keV_erg**3 / h_plank_cgs**3 / c_light_cgs**2

    exp_arg = energy_keV / (f_color_factor * kbTeff_keV)
    # exp_arg = np.float128(exp_arg)
    # # Cap the argument to np.exp() to avoid overflow
    # exp_arg = np.minimum(
    #     exp_arg, 1000
    # )  # np.exp(700) is roughly the upper limit before overflow
    warnings.filterwarnings("ignore")
    return A1 * energy_keV**3 * Y_angle_func / (np.exp(exp_arg) - 1.0)


def FluxDimless_To_kbTeff_keV(F_dimless, Mass_Msun, Mdot):
    Msun = 1
    A2 = 0.1331 * np.power((1.0e18 / Mdot), 0.25) * np.power(Mass_Msun / Msun, 0.5)

    return np.power(F_dimless, 0.25) / A2


def FindIsco(a):
    # accounts for negative spin
    sign = 1.0
    if a < 0.0:
        sign = -1.0

    Z1 = 1.0 + math.pow(1.0 - a * a, 1.0 / 3.0) * (
        math.pow(1.0 + a, 1.0 / 3.0) + math.pow(1.0 - a, 1.0 / 3.0)
    )
    Z2 = math.sqrt((3.0 * a * a) + (Z1 * Z1))

    return 3.0 + Z2 - sign * math.sqrt((3.0 - Z1) * (3.0 + Z1 + (2 * Z2)))


def KerrEcirc(r, a):
    u = 1 / r
    u2 = u * u
    u3 = u2 * u

    A = 1.0 - 2.0 * u + a * math.sqrt(u3)
    B = math.sqrt(1.0 - 3.0 * u + 2.0 * a * math.sqrt(u3))

    return A / B


def glpRed(r, h, a):

    sqtr = r**0.5
    r2 = r * r
    h2 = h * h
    a2 = a * a
    return (
        (r * sqtr + a)
        * (h2 - 2 * h + a2) ** 0.5
        / sqtr
        / (r2 - 3 * r + 2 * a * sqtr) ** 0.5
        / (h2 + a2) ** 0.5
    )


def find_isco(a):
    # accounts for negative spin
    sign = 1.0
    if a < 0.0:
        sign = -1.0
    Z1 = 1.0 + math.pow(1.0 - a * a, 1.0 / 3.0) * (
        math.pow(1.0 + a, 1.0 / 3.0) + math.pow(1.0 - a, 1.0 / 3.0)
    )
    Z2 = math.sqrt((3.0 * a * a) + (Z1 * Z1))
    return 3.0 + Z2 - sign * math.sqrt((3.0 - Z1) * (3.0 + Z1 + (2 * Z2)))


def Lorens(r, a):
    r2 = r * r
    a2 = a * a
    r3 = r2 * r
    r1_4 = r**0.25
    r3_2 = r3**0.5
    r1_2 = r**0.5

    return (
        (r2 - 2 * r + a2) ** 0.5
        * (r3_2 + a)
        / r1_4
        / (r * r1_2 + 2 * a - 3 * r1_2) ** 0.5
        / (r3 + a2 * r + 2 * a2) ** 0.5
    )


def dA_dr(r, a):

    r2 = r * r
    a2 = a * a
    r3 = r2 * r
    r4 = r3 * r

    # d1 = r4 + a2*r2 + 2*a2*r
    d1 = r2 + a2 + 2 * a2 / r
    d2 = r2 - 2 * r + a2

    return 2 * np.pi * r * (d1 / d2) ** 0.5


# Novikov-Thorne Disk
###############################################

chi = None


def gmetric(r, th):
    m = np.zeros((4, 4))
    t3 = chi**2
    t2 = r**2
    t5 = np.cos(th)
    t6 = t5**2
    t7 = t3 * t6
    t8 = t2 + t7
    t9 = 1 / t8
    t19 = np.sin(th)
    t20 = t19**2

    gtt = -1 + 2 * r * t9
    grr = t8 / ((-2 + r) * r + t3)
    gthth = t8
    gpp = t20 * (t2 + t3 + 2 * r * t20 * t3 * t9)
    gtp = -2 * chi * r * t20 * t9

    m[0][0] = gtt
    m[0][3] = gtp
    m[1][1] = grr
    m[2][2] = gthth
    m[3][0] = gtp
    m[3][3] = gpp

    return m


def metric_dr_dr2(r, th):
    m = gmetric(r, th)
    dg = np.zeros((4, 4))
    dg2 = np.zeros((4, 4))

    t3 = chi**2
    t2 = r**2
    t5 = np.cos(th)
    t6 = t5**2
    t7 = t3 * t6
    t8 = t2 + t7
    t9 = 1 / t8
    t19 = np.sin(th)
    t20 = t19**2
    t25 = t8**-2
    t38 = r**3
    t39 = t8**-3

    dgttdr = -4 * t2 * t25 + 2 * t9
    dgppdr = t20 * (2 * r - 4 * t2 * t20 * t25 * t3 + 2 * t20 * t3 * t9)
    dgtpdr = 4 * chi * t2 * t20 * t25 - 2 * chi * t20 * t9
    dgttdr2 = -12 * r * t25 + 16 * t38 * t39
    dgppdr2 = t20 * (2 - 12 * r * t20 * t25 * t3 + 16 * t20 * t3 * t38 * t39)
    dgtpdr2 = 12 * chi * r * t20 * t25 - 16 * chi * t20 * t38 * t39

    dg[0][0] = dgttdr
    dg[0][3] = dgtpdr
    dg[3][0] = dgtpdr
    dg[3][3] = dgppdr

    dg2[0][0] = dgttdr2
    dg2[0][3] = dgtpdr2
    dg2[3][0] = dgtpdr2
    dg2[3][3] = dgppdr2

    return m, dg, dg2


def IntegFunc(r, alpha):
    m, dmdr, dmdr2 = metric_dr_dr2(r, pi / 2)

    Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3] ** 2 - dmdr[0][0] * dmdr[3][3])) / dmdr[
        3
    ][3]
    denom = sqrt(-(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var**2))
    E_var = -(m[0][0] + m[0][3] * Omega_var) / denom
    Lz_var = (m[0][3] + m[3][3] * Omega_var) / denom

    dOmega_var_dr = (
        -Omega_var * dmdr2[3][3] / dmdr[3][3]
        + (
            -dmdr2[0][3]
            + (
                2 * dmdr[0][3] * dmdr2[0][3]
                - dmdr2[0][0] * dmdr[3][3]
                - dmdr[0][0] * dmdr2[3][3]
            )
            / 2
            / sqrt(dmdr[0][3] ** 2 - dmdr[0][0] * dmdr[3][3])
        )
        / dmdr[3][3]
    )

    exp1 = (
        -dmdr[0][0]
        - 2 * (dmdr[0][3] * Omega_var + m[0][3] * dOmega_var_dr)
        - (dmdr[3][3] * Omega_var**2 + m[3][3] * 2 * Omega_var * dOmega_var_dr)
    )
    exp2 = dmdr[0][3] + dmdr[3][3] * Omega_var + m[3][3] * dOmega_var_dr

    dLz_var_dr = -Lz_var / (2 * denom**2) * exp1 + exp2 / denom

    return (E_var - Omega_var * Lz_var) * dLz_var_dr


def IntegFlux(re, rstart):
    alpha = 1.0
    result, error = quad(IntegFunc, rstart, re, args=(alpha,))
    return result


def FluxDimless(spin, re, rstart):
    global chi
    chi = spin
    m, dmdr, dmdr2 = metric_dr_dr2(re, pi / 2)

    Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3] ** 2 - dmdr[0][0] * dmdr[3][3])) / dmdr[
        3
    ][3]
    denom = sqrt(-(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var**2))
    E_var = -(m[0][0] + m[0][3] * Omega_var) / denom
    Lz_var = (m[0][3] + m[3][3] * Omega_var) / denom

    dOmega_var_dr = (
        -Omega_var * dmdr2[3][3] / dmdr[3][3]
        + (
            -dmdr2[0][3]
            + (
                2 * dmdr[0][3] * dmdr2[0][3]
                - dmdr2[0][0] * dmdr[3][3]
                - dmdr[0][0] * dmdr2[3][3]
            )
            / 2
            / sqrt(dmdr[0][3] ** 2 - dmdr[0][0] * dmdr[3][3])
        )
        / dmdr[3][3]
    )

    integ = IntegFlux(re, rstart)
    sqdetg = sqrt(m[1][1] * (m[0][3] ** 2 - m[0][0] * m[3][3]))

    return -1.0 / sqdetg * dOmega_var_dr / (E_var - Omega_var * Lz_var) ** 2 * integ


######################  End Novikov-Thorne disk


def SaveXY(X, Y, fname):
    with open(fname, "w") as f:
        for x, y in zip(X, Y):
            txt = "%E %E\n" % (x, y)
            f.write(txt)


def FluxFromSpectra(spectra, energy):
    erg = 1
    keV = erg / 624150647.99632
    Flux = cumtrapz(spectra, energy, initial=0)[-1] * keV
    return Flux


def IonazationFromFlux(Flux, density):
    return 4 * np.pi * Flux / density



def ln_power_scale_r_array_1_2(N, rmin, rmax):
    """
    Generate an array with spacing according to (ln(x))^1.2.

    Parameters:
    N (int): Number of points
    rmin (float): Minimum value (start of the array)
    rmax (float): Maximum value (end of the array)

    Returns:
    numpy.ndarray: Array of values spaced according to (ln(x))^1.2.
    """
    # Generate N points linearly spaced between ln(rmin) and ln(rmax)
    ln_space = np.linspace(np.log(rmin), np.log(rmax), N)
    
    # Apply the transformation (ln(x))^1.2
    r_array = np.exp(ln_space**1.2)
    
    # Normalize to ensure the values are within the range [rmin, rmax]
    r_array = (r_array - r_array.min()) / (r_array.max() - r_array.min()) * (rmax - rmin) + rmin
    
    return r_array


def GenGlobalRadiArray(spin, N):
    rmin = FindIsco(spin) + 1.0e-5
    rmax = 500.
    Nradi = N
    r_array = ln_power_scale_r_array_1_2(Nradi, rmin, rmax)
    r_array_reversed = r_array[::-1]
    filename = "data/data_local/global_data/radi_a%.3f.txt" % spin
    np.savetxt(filename, r_array_reversed)

    return r_array_reversed