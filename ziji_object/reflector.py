import os
import numpy as np
from scipy.integrate import cumtrapz


from .parallel_executor import ParallelExecutor

from .ziji_param import AccretionDisk
from .ziji_utils import SaveXY


class Reflector:
    def __init__(self, accretion_disk: AccretionDisk, mode=0):
        # def __init__(self, ne_density_arr, ionazation_arr, A_fe, spec_irradiators, mode=0):
        # if mode 0 n_density_arr and ionazation_arr are just number
        self.ne_density_arr = accretion_disk.density
        self.ionazation_arr = accretion_disk.ionization  # Xi_arr
        self.mode = mode  # 0 - constant Xi, 1 - self consistent Xi
        self._init_profile(accretion_disk)
        self.A_fe = accretion_disk.iron_abundance

        self.spec_irradiators = None
        ##################

        self.spec_reflections = None
        self.ionization_min = 1.0
        self.emis_prof_eff = None

    def _init_profile(self, accretion_disk: AccretionDisk):
        if self.mode == 1:
            self.ne_density_arr = accretion_disk.density_profile
            # self.ionazation_arr = accretion_disk.ionization_profile

    def SaveSpec(self, iter=0):
        fname = (
            "data/data_out/spectra_accretion_disk/"
            + "iter%d/reflection_iter%d.npy" % (iter, iter)
        )
        np.save(fname, self.spec_reflections)

        fname = (
            "data/data_out/spectra_accretion_disk/"
            + "iter%d/emis_prof%d.txt" % (iter, iter)
        )
        SaveXY(np.arange(100), self.emis_prof_eff, fname)

    def ParallelReflionx(self, spec_irradiators: np.ndarray, num_processes=100):
        self.spec_irradiators = spec_irradiators
        self.spec_reflections = np.zeros(spec_irradiators.shape)

        n_irradiators, n_energy_arr = (self.spec_irradiators).shape
        energy_reflionx = np.loadtxt(
            "data/data_local/global_data/energybins", unpack=True
        )

        ##INIT reflionx
        self._CreatInputSpec(energy_reflionx, n_irradiators)
        self.emis_prof_eff = self._EmisProfArr(energy_reflionx)
        self._CreatInputScript(n_irradiators)

        # run reflionx
        binary_commands = self._BinaryCommands(n_irradiators)
        executor_binary = ParallelExecutor(
            num_processes, mode="binary", commands=binary_commands
        )
        executor_binary.execute()

        # Packing data
        self._PackOutSpec(n_irradiators)
        self.DeleteCache()

        return self.spec_reflections

    def _BinaryCommands(self, n_commands):
        list_commands = [
            [
                "lib/reflector_linux",
                f"data/data_cache/reflionx/input_for_script/r{i}",
            ]
            for i in range(n_commands)
        ]
        return list_commands

    def _EmisProfArr(self, energy_inc_arr):
        n_irradiators, n_energy_arr = (self.spec_irradiators).shape
        FxArr_eff = np.zeros(n_irradiators)
        erg = 1
        keV = erg / 624150647.99632
        for i in range(n_irradiators):
            FxArr_eff[i] = (
                cumtrapz(self.spec_irradiators[i, :], energy_inc_arr, initial=0)[-1]
                * keV
            )

        Xi_prof_eff = 4 * np.pi * FxArr_eff / self.ne_density_arr

        if self.mode == 0:
            emis_eff = Xi_prof_eff / self.ionazation_arr
            return emis_eff
        else:
            emis_eff = np.full(100, 1.)
            self.ionazation_arr = Xi_prof_eff
            boolmask = Xi_prof_eff < self.ionization_min
            # if np.all(boolmask):
            emis_eff[boolmask] = Xi_prof_eff[boolmask] / self.ionization_min

            return emis_eff

    def _CreatInputSpec(self, energy_arr, n_irradiators):
        for i in range(n_irradiators):
            fname = "data/data_cache/reflionx/initial_spec/r%d.txt" % i
            SaveXY(energy_arr, self.spec_irradiators[i, :], fname)

    def _PackOutSpec(self, n_irradiators):
        for i in range(n_irradiators):
            En, c_E_F_E = np.loadtxt(
                "data/data_cache/reflionx/out_spec/r%d.txt" % i, unpack=True
            )
            # c_F_E = c_E_F_E / En
            # erg = 1
            # keV = erg / 624150647.99632
            # Fx = cumtrapz(c_F_E, En, initial=0)[-1] * keV
            # if self.mode == 1:
            #     Xi = 4 * np.pi * Fx / self.ne_density_arr[i]
            #     self.spec_reflections[i, :] = c_F_E * self.ionazation_arr[i]/Xi
            # elif self.mode == 0:
            #     Xi = 4 * np.pi * Fx / self.ne_density_arr
            #     self.spec_reflections[i, :] = c_F_E * self.ionazation_arr/Xi

            self.spec_reflections[i, :] = self.emis_prof_eff[i] * 10**20 * c_E_F_E / En

    def _CreatInputScript(self, n_irradiators):
        Xi = self.ionazation_arr
        ne_dens = self.ne_density_arr
        for i in range(n_irradiators):
            fname = "data/data_cache/reflionx/input_for_script/r%d" % i
            if self.mode == 1:
                Xi = self.ionazation_arr[i]
                ne_dens = self.ne_density_arr[i]
                if Xi < self.ionization_min:
                    Xi = self.ionization_min

            with open(fname, "w") as f:
                f.write("1\n")
                f.write("1 'data/data_cache/reflionx/out_spec/r%d.txt'\n" % i)
                f.write(
                    "%E %E %.1f 'data/data_cache/reflionx/initial_spec/r%d.txt'\n"
                    % (ne_dens, Xi, self.A_fe, i)
                )
                f.write(
                    "'data/data_cache/reflionx/other_out_data/refbhb_0_output_r%d.out'\n"
                    % i
                )
                f.write(
                    "'data/data_cache/reflionx/other_out_data/refbhb_0_output_r%d.dat'\n"
                    % i
                )
                # f.write(
                #     "'data/data_cache/reflionx/other_out_data/refbhb_0_spectra_r%d.dat'\n"
                #     % i
                # )

    def DeleteCache(self):
        directory = "data/data_cache/reflionx/"
        # Walk through the directory
        for root, dirs, files in os.walk(directory):
            for file in files:
                file_path = os.path.join(root, file)
                os.remove(file_path)  # Delete the file
