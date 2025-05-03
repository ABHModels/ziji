import numpy as np
import mpmath
from scipy.integrate import cumtrapz
import math


from .ziji_param import Corona, BlackHole, SpectralDiskGrid
from .ziji_utils import (
    glpRed,
)
from .ziji_utils import SaveXY


class Lampost:
    # __init__(self, height, spin, gamma, ecut, lp_dat)
    def __init__(self, corona: Corona, black_hole: BlackHole):
        self.height = corona.height
        self.gamma = corona.gamma
        self.ecut = corona.ecut

        self.scale_luminosity = corona.scale_luminosity
        ###################

        # black hole
        self.mass_black_hole = black_hole.mass
        self.spin = black_hole.spin

        self.lp_dat = None  # lp_dat
        self._load_data()

        self.elow = 0.001 #1eV
        self.ehigh = mpmath.inf

        self.coeff_4pi = None  # Adam ingram, reltrans
        self.Fx_inner = None  # in erg
        # self.radi_irradiation = None
        # self.energy_arr = None
        self.emis_profile = None  # dimless
        self.flux_profile = None  # in erg
        self.spec_irradiation = None  # 2d numpy array, F_E
        self.direct_spec = None

    def _load_data(self):
        self.lp_dat = np.loadtxt(
            "data/data_local/lampost/lp_%.3f_%.3f.dat" % (self.spin, self.height),
            unpack=True,
        )

    def SaveSpec(self):
        np.save(
            "data/data_out/spectra_accretion_disk/lampost_spec.npy",
            self.spec_irradiation,
        )

    # energy_arr in keV
    def NormWithXi(self, Xi_in, ndens_in, spectral_grid_irradiation: SpectralDiskGrid):

        radi_arr = spectral_grid_irradiation.radii
        energy_arr = spectral_grid_irradiation.energies

        index_rmin = np.argmin(radi_arr)
        radipr, Inhr = self.lp_dat
        Inhr_arr = np.interp(radi_arr, radipr, Inhr)
        glp_arr = glpRed(radi_arr, self.height, self.spin)
        glpInten_arr = glp_arr**self.gamma * (2.0 * np.pi * Inhr_arr)
        glpInten_arr_integ = glp_arr**2 * (2.0 * np.pi * Inhr_arr)

        elow, ehigh = np.min(energy_arr), np.max(energy_arr)
        keV = 1
        erg = 624150647.99632 * keV
        # LpSp_Ximax = self._SpecFlux_CutoffPow( energy_arr, self.ecut*glp_arr[index_rmin])
        # IntegPow = glpInten_arr[index_rmin]*cumtrapz(LpSp_Ximax, energy_arr, initial=0)[-1]
        Integ_Flux = self._IntegratedFlux(self.gamma, self.ecut, self.elow, self.ehigh)
        IntegPow = glpInten_arr_integ[index_rmin] * Integ_Flux
        self.Fx_inner = Xi_in * ndens_in / (4 * np.pi) * erg  # Fx in keV now
        self.coeff_4pi = self.Fx_inner / IntegPow

        Nradi, Nspec = len(radi_arr), len(energy_arr)
        # lp_spec = np.zeros((Nradi, Nspec))
        emis_prof = np.zeros(Nradi)

        for i in range(Nradi):
            spectral_grid_irradiation.spectra[i, :] = (
                self.coeff_4pi
                * glpInten_arr[i]
                * self._SpecFlux_CutoffPow(energy_arr, self.ecut * glp_arr[i])
            )
            Integ_Flux_emis = self._IntegratedFlux(
                self.gamma, glp_arr[i] * self.ecut, elow, ehigh
            )
            emis_prof[i] = (
                1.0 / self.Fx_inner * self.coeff_4pi * glpInten_arr[i] * Integ_Flux_emis
            )

        self.radi_irradiation = radi_arr
        # self.energy_arr = energy_arr
        self.emis_profile = emis_prof
        self.flux_profile = (
            self.coeff_4pi * glpInten_arr_integ * Integ_Flux / erg
        )  # Flux in erg now

        self.spec_irradiation = spectral_grid_irradiation.spectra
        self.Fx_inner = self.Fx_inner / erg  # Flux in erg now

        return spectral_grid_irradiation

    # energy_arr in keV
    # Luminosity in erg/s
    # corona_luminosity_M2 ----- Luminosity/(GM/c^2)^2
    def NormWithCoronaLuminosity(self, spectral_grid_irradiation: SpectralDiskGrid):
        ratio_lum_edd = self.scale_luminosity
        Mass_msun = self.mass_black_hole

        radi_arr = spectral_grid_irradiation.radii
        energy_arr = spectral_grid_irradiation.energies

        index_rmin = np.argmin(radi_arr)
        radipr, Inhr = self.lp_dat
        Inhr_arr = np.interp(radi_arr, radipr, Inhr)
        glp_arr = glpRed(radi_arr, self.height, self.spin)
        glpInten_arr = glp_arr**self.gamma * (2.0 * np.pi * Inhr_arr)
        glpInten_arr_integ = glp_arr**2 * (2.0 * np.pi * Inhr_arr)

        elow, ehigh = np.min(energy_arr), np.max(energy_arr)

        Integ_Flux = self._IntegratedFlux(self.gamma, self.ecut, self.elow, self.ehigh)
        IntegPow = glpInten_arr_integ[index_rmin] * Integ_Flux

        keV = 1
        erg = 624150647.99632 * keV
        G_CGS = 6.67428e-8  # Gravitational constant (cgs: cm^3 * g^-1 * s-2)
        MASS_SUN_CGS = 1.98843e33  # g
        C_LIGHT = 29979245800.0  # cm/s
        r_grav = G_CGS * Mass_msun * MASS_SUN_CGS / C_LIGHT**2
        Led_M = 1.26e38  # erg/s  1/Msun
        luminosity = ratio_lum_edd * Led_M * Mass_msun  # erg/s
        corona_luminosity_M2 = luminosity / r_grav**2

        self.coeff_4pi = (corona_luminosity_M2 * erg) / (4 * np.pi) / Integ_Flux
        self.Fx_inner = self.coeff_4pi * IntegPow

        Nradi, Nspec = len(radi_arr), len(energy_arr)
        # lp_spec = np.zeros((Nradi, Nspec))
        emis_prof = np.zeros(Nradi)

        for i in range(Nradi):
            spectral_grid_irradiation.spectra[i, :] = (
                self.coeff_4pi
                * glpInten_arr[i]
                * self._SpecFlux_CutoffPow(energy_arr, self.ecut * glp_arr[i])
            )
            Integ_Flux_emis = self._IntegratedFlux(
                self.gamma, glp_arr[i] * self.ecut, elow, ehigh
            )
            emis_prof[i] = (
                1.0 / self.Fx_inner * self.coeff_4pi * glpInten_arr[i] * Integ_Flux_emis
            )

        self.radi_irradiation = radi_arr
        # self.energy_arr = energy_arr
        self.emis_profile = emis_prof
        self.flux_profile = (
            self.coeff_4pi * glpInten_arr_integ * Integ_Flux / erg
        )  # Flux in erg now
        self.spec_irradiation = spectral_grid_irradiation.spectra
        self.Fx_inner = self.Fx_inner / erg  # Flux in erg now

        return spectral_grid_irradiation

    # dOmegaDet_dADet -- dimless, direct spec -- F_E
    def CalculateDirectSpec(self, inc, distance, energy_arr):
        dOmegaDet_dADet, redshift_source_observer = self._InitDirec(
            int(inc), int(self.height)
        )
        dOmegaDet_dADet = dOmegaDet_dADet / distance**2
        f_E = self._SpecFlux_CutoffPow(energy_arr, self.ecut * redshift_source_observer)
        self.direct_spec = (
            self.coeff_4pi
            * redshift_source_observer**self.gamma
            * f_E
            * dOmegaDet_dADet
        )

        fname = f"data/data_out/spectra_observed/corona_inc{inc:.0f}.txt"
        SaveXY(energy_arr, self.direct_spec, fname)

    def _InitDirec(self, inc_int, h_int):
        direc_dic = {}
        direc_dic[45] = {
            2: [1.999263e-01, 4.464976e-01],
            3: [4.004033e-01, 6.322659e-01],
            5: [6.1607e-01, 7.844269e-01],
            10: [8.032126e-01, 8.955291e-01],
            20: [9.028092e-01, 9.488142e-01]
        }

        return direc_dic[inc_int][h_int]

    # F_E
    def _SpecFlux_CutoffPow(self, E, gecut):
        return np.power(E, (1.0 - self.gamma)) * np.exp(-E / gecut)

    # in keV = 1
    def _IntegratedFlux(self, gamma, ecut, elow, ehigh):
        Flux = np.power(ecut, (2 - gamma)) * mpmath.gammainc(
            2.0 - gamma, elow / ecut, ehigh / ecut
        )
        return float(Flux)  # because of mpmath format

    def _EmisProfArr(self, Xi_in, radi_arr, ndens_in, energy_arr, irradi_spec2d):

        Nradi = len(radi_arr)
        FxArr = np.zeros(Nradi)
        erg = 1
        keV = erg / 624150647.99632
        for i in range(Nradi):
            FxArr[i] = cumtrapz(irradi_spec2d[i, :], energy_arr, initial=0)[-1] * keV

        XiArr = 4 * np.pi * FxArr / ndens_in

        return XiArr / Xi_in
