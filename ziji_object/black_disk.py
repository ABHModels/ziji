import numpy as np

from .ziji_param import BlackHole, AccretionDisk, SpectralDiskGrid

from .ziji_utils import (
    BlackBodyIntensity,
    FluxDimless_To_kbTeff_keV,
    FindIsco,
    KerrEcirc,
    FluxDimless,
)

# Novikov-Thorne Disk


class BlacDisk:
    # __init__(self, ratio_lum_edd, Mass_msun, spin, f_color_factor)
    def __init__(self, black_hole: BlackHole, accretion_disk: AccretionDisk):
        self.Mass_msun = black_hole.mass  # in Msun, dimless
        self.ratio_lum_edd = (
            accretion_disk.scale_luminosity  # in L_edington, dimless, for example 0.1 * Led
        )
        self.spin = black_hole.spin

        self.Y_angle_func = 1.0  # isotropic emission
        self.f_color_factor = accretion_disk.color_factor
        ##########################

        self.luminosity = None
        self.mass_acc_rate = None
        self._init_addit_param()

        # self.radi_emitters = None  # 1d array
        # self.energy_arr = None
        self.spec_emitters = None  # 2d array, F_E

    def SaveSpec(self):
        np.save(
            "data/data_out/spectra_accretion_disk/black_disk_spec.npy",
            self.spec_emitters,
        )

    # energy_arr in keV
    def InitBlackSpec(self, spectral_grid_emitters: SpectralDiskGrid):
        radi_arr = spectral_grid_emitters.radii
        energy_arr = spectral_grid_emitters.energies
        Nradi, Nspec = len(radi_arr), len(energy_arr)
        # bb_spec = np.zeros((Nradi, Nspec))

        isco = FindIsco(self.spin)
        for i in range(Nradi):
            F_dimless = FluxDimless(self.spin, radi_arr[i], isco)
            kbTeff_keV = FluxDimless_To_kbTeff_keV(
                F_dimless, self.Mass_msun, self.mass_acc_rate
            )
            spectral_grid_emitters.spectra[i, :] = np.pi * BlackBodyIntensity(
                energy_arr, kbTeff_keV, self.f_color_factor, self.Y_angle_func
            )

        self.spec_emitters = spectral_grid_emitters.spectra
        return spectral_grid_emitters
        # self.radi_emitters = radi_arr
        # self.energy_arr = energy_arr
        # self.spec_emitters = bb_spec

    def _init_addit_param(self):
        Led_M = 1.26e38  # erg/s  1/Msun
        self.luminosity = self.ratio_lum_edd * Led_M * self.Mass_msun
        isco = FindIsco(self.spin)
        E_isco = KerrEcirc(isco, self.spin)
        epsil = 1.0 - E_isco
        C_LIGHT = 29979245800.0  # cm/s
        self.mass_acc_rate = self.luminosity / (epsil * C_LIGHT**2)
