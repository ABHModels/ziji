import numpy as np


from . import (
    LampostCorona,
    BlacDisk,
    ConvolReturn,
    Reflector,
    SpectralTransit,
    ZijiParams,
    SpectralDiskGrid,
)

from .ziji_utils import FluxFromSpectra, IonazationFromFlux


class Iterator:
    def __init__(self, ziji_params: ZijiParams, n_iter=4, mode=0):
        self.ziji_params = ziji_params
        self.n_iter = n_iter
        self.mode = mode

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
        lp_spec_grid = lampost_corona.NormWithCoronaLuminosity(lp_spec_grid)
        lampost_corona.SaveSpec()
        self.corona_irraditon = lp_spec_grid.spectra
        # lampost_corona.CalculateDirectSpec(
        #     self.ziji_params.observer.inclination,
        #     self.ziji_params.observer.distance,
        #     self.ENERGYS,
        # )

        black_disk_grid = SpectralDiskGrid(radii=self.RADII, energies=self.ENERGYS)
        black_disk = BlacDisk(
            self.ziji_params.black_hole, self.ziji_params.accretion_disk
        )
        black_disk_grid = black_disk.InitBlackSpec(black_disk_grid)
        black_disk.SaveSpec()
        self.black_disk_emission = black_disk_grid.spectra

        spec_transit = SpectralTransit(
            self.ziji_params.black_hole, self.ziji_params.observer
        )
        spec_transit.Transit(black_disk_grid, self.ENERGYS)
        inc = self.ziji_params.observer.inclination
        spec_transit.SaveSpec(f"black_disk_inc{inc:.0f}.txt")

        return_black_disk_grid = SpectralDiskGrid(
            radii=self.RADII, energies=self.ENERGYS
        )
        convol_return = ConvolReturn(
            self.ziji_params.black_hole,
            return_black_disk_grid,
            black_disk_grid,
        )
        return_black_disk_grid = convol_return.ParallelConvol(num_processes=100)
        convol_return.SaveSpec()
        self.returning_black = return_black_disk_grid.spectra

        flux_return_inner = FluxFromSpectra(self.returning_black[-1, :], self.ENERGYS)
        flux_inner = lampost_corona.Fx_inner + flux_return_inner
        Xi_inner = IonazationFromFlux(
            flux_inner, self.ziji_params.accretion_disk.density_profile[-1]
        )
        self.ziji_params.accretion_disk.ionization = Xi_inner

        return Xi_inner

    def Loop(self):
        for i in range(self.n_iter + 1):
            self._Iter(i)

    def _Iter(self, iter):
        if iter == 0:
            reflector = Reflector(self.ziji_params.accretion_disk, mode=self.mode)
            spec_reflections = reflector.ParallelReflionx(
                self.corona_irraditon + self.returning_black
            )
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

            total_irradiation = (
                self.corona_irraditon
                + self.returning_black
                + return_reflection_grid.spectra
            )
            reflector = Reflector(self.ziji_params.accretion_disk, mode=self.mode)
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
