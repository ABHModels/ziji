import numpy as np
from scipy.interpolate import RegularGridInterpolator as reginter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor


from .ziji_param import BlackHole, SpectralDiskGrid


class ConvolReturn:
    def __init__(
        self,
        black_hole: BlackHole,
        spectral_grid_refer: SpectralDiskGrid,
        spectral_grid_emitters: SpectralDiskGrid,
    ):
        # def __init__(
        #     self,
        #     spin,
        #     radi_refer_points,
        #     energy_refer_points,
        #     radi_emitters,
        #     energy_emitters,
        #     spec_emitters,
        # ):

        self.spin = black_hole.spin
        self.radi_refer_points = spectral_grid_refer.radii
        self.energy_refer_points = spectral_grid_refer.energies

        self.radi_emitters = spectral_grid_emitters.radii
        self.energy_emitters = spectral_grid_emitters.energies
        self.spec_emitters = spectral_grid_emitters.spectra
        ###########################

        self.n_refer_points = len(self.radi_refer_points)
        self.spec_refer_points = np.zeros(
            (self.n_refer_points, len(self.energy_refer_points))
        )

    def ParallelConvol(self, num_processes=100):
        k_refer_points_list = [i for i in range(self.n_refer_points)]
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            # Collect results
            results = executor.map(self._GenRetDat, k_refer_points_list)
            # Update the main object's state based on the results
            for k, result in enumerate(results):
                self.spec_refer_points[k, :] = result

        return SpectralDiskGrid(
            radii=self.radi_refer_points,
            energies=self.energy_refer_points,
            spectra=self.spec_refer_points,
        )

    def SaveSpec(self, iter=None):
        if iter is None:
            fname = "data/data_out/spectra_accretion_disk/return_black_disk.npy"
        else:
            fname = (
                "data/data_out/spectra_accretion_disk/"
                + "iter%d/return_reflection_iter%d.npy" % (iter, iter)
            )

        np.save(fname, self.spec_refer_points)

    def _GenRetDat(self, k_refer_point):

        _, dphiArr, dlambArr, _, _, _, _ = np.loadtxt(
            "data/data_local/returning_radiation/robInPar_a%.3f.txt" % (self.spin),
            unpack=True,
        )
        fnrrad = "data/data_local/returning_radiation/radrem_a%.3f/r%d.npy" % (
            self.spin,
            k_refer_point,
        )
        radrem = np.load(fnrrad)
        dldph = dphiArr[k_refer_point] * dlambArr[k_refer_point]
        RetRadSpec_ir = dldph * self._SumIobsRetRadSpec(radrem)
        # np.save('data/r%d.npy'%k_refer_point, RetRadSpec_ir)
        return RetRadSpec_ir

    def _DelSort(self, radrem):
        rarr1d = radrem[:, :, 0].flatten()
        garr1d = radrem[:, :, 1].flatten()
        # N = len(rarr1d)

        rarrSort = rarr1d[rarr1d > 1.0e-10]
        garrSort = garr1d[rarr1d > 1.0e-10]

        return rarrSort, garrSort

    def _FemArr_InterpolSpec(self, Eob, rarr, garr, Eem, radiSpec, RefSpec):

        EemArr = Eob / garr
        # SpeArr = np.interp(EemArr, Eem, Spem)
        fin = reginter((radiSpec, Eem), RefSpec, bounds_error=False, fill_value=0)
        FemArr = fin((rarr, EemArr))
        # SpeArr[np.logical_or(EemArr<Emin, EemArr>Emax)] = 0.
        return FemArr  # [np.logical_xor(EemArr<Emin, EemArr>Emax)]

    def _SumIobsRetRadSpec(self, radrem_retrad):  # Spem -> F_E
        rarrS, garrS = self._DelSort(radrem_retrad)
        glp3_pi = garrS**3 / np.pi  # Iobs = g^3 Iem = g^3 Fem/pi

        Neobs = len(self.energy_refer_points)
        SpecRetRad = np.zeros(Neobs)
        for i, Eob in enumerate(self.energy_refer_points):
            FemArr = self._FemArr_InterpolSpec(
                Eob,
                rarrS,
                garrS,
                self.energy_emitters,
                self.radi_emitters,
                self.spec_emitters,
            )
            IobsRetRadArr = glp3_pi * FemArr
            SpecRetRad[i] = np.sum(IobsRetRadArr)

        return SpecRetRad  # Sum Specific Intensity for all phot
