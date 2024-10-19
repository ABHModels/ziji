import numpy as np


from . import (
    LampostCorona,
    ConvolReturn,
    Reflector,
    SpectralTransit,
    ZijiParams,
    SpectralDiskGrid,
)


class Iterator:
    def __init__(self, ziji_params: ZijiParams, n_iter=4):
        self.ziji_params = ziji_params
        self.n_iter = n_iter

        self.ENERGYS = np.loadtxt("data/data_local/global_data/energybins", unpack=True)
        self.RADII = np.loadtxt(
            "data/data_local/global_data/radi_a%.3f.txt" % ziji_params.black_hole.spin,
            unpack=True,
        )

        self.tmp_reflection_grid = None
        self.tmp_returning = None  # Total
        self.returning_black = None  # Black Disk returning
        self.black_disk_emission = None
        self.corona_irraditon = None

    def InitModel(self):
        lp_spec_grid = SpectralDiskGrid(radii=self.RADII, energies=self.ENERGYS)
        lampost_corona = LampostCorona(
            self.ziji_params.corona, self.ziji_params.black_hole
        )
        # lp_spec_grid = lampost_corona.NormWithCoronaLuminosity(lp_spec_grid)
        lp_spec_grid = lampost_corona.NormWithXi(
            self.ziji_params.accretion_disk.ionization,
            self.ziji_params.accretion_disk.density,
            lp_spec_grid,
        )
        lampost_corona.SaveSpec()

        self.corona_irraditon = lp_spec_grid.spectra

    def Loop(self):
        for i in range(self.n_iter + 1):
            self._Iter(i)

    def _Iter(self, iter):
        if iter == 0:
            reflector = Reflector(self.ziji_params.accretion_disk)
            spec_reflections_lp = reflector.ParallelReflionx(self.corona_irraditon)
            reflector.SaveSpec(iter=iter)
            self.tmp_reflection_grid = SpectralDiskGrid(
                radii=self.RADII,
                energies=self.ENERGYS,
                spectra=spec_reflections_lp,
            )

            inc = self.ziji_params.observer.inclination
            spec_transit = SpectralTransit(
                self.ziji_params.black_hole, self.ziji_params.observer
            )
            spec_transit.Transit(self.tmp_reflection_grid, self.ENERGYS)
            spec_transit.SaveSpec(f"reflection_inc{inc:.0f}_iter{iter}.txt")

        else:
            return_reflection_grid = SpectralDiskGrid(
                radii=self.RADII, energies=self.ENERGYS
            )

            convol_return = ConvolReturn(
                self.ziji_params.black_hole,
                return_reflection_grid,
                self.tmp_reflection_grid,
            )
            return_reflection_grid = convol_return.ParallelConvol(num_processes=100)
            convol_return.SaveSpec(iter=iter)

            total_irradiation = self.corona_irraditon + return_reflection_grid.spectra
            reflector = Reflector(self.ziji_params.accretion_disk)
            spec_reflections = reflector.ParallelReflionx(total_irradiation)
            reflector.SaveSpec(iter=iter)
            self.tmp_reflection_grid = SpectralDiskGrid(
                radii=self.RADII,
                energies=self.ENERGYS,
                spectra=spec_reflections,
            )

            inc = self.ziji_params.observer.inclination
            spec_transit = SpectralTransit(
                self.ziji_params.black_hole, self.ziji_params.observer
            )
            spec_transit.Transit(self.tmp_reflection_grid, self.ENERGYS)
            spec_transit.SaveSpec(f"reflection_inc{inc:.0f}_iter{iter}.txt")
