import numpy as np
from scipy.interpolate import RegularGridInterpolator as reginter


from .ziji_param import BlackHole, Observer, SpectralDiskGrid
from .ziji_utils import SaveXY


# The spectra moving or transitioning from the disk to the observer
class SpectralTransit:
    def __init__(
        self,
        black_hole: BlackHole,
        observer: Observer,
        # transit_data_dict,
    ):
        # def __init__(
        #     self,
        #     spin,
        #     inc_angle,
        #     distance,
        #     transit_data_dict,
        # ):
        self.spin = black_hole.spin
        self.inc_angle = observer.inclination
        self.distance = observer.distance  # in Mass, dimless

        self.transit_data_dict = None  # transit_data_dict
        self._load_data()
        ###############################

        self.radi_disk_emitters = None
        self.energy_disk_emitters = None
        self.spec_disk_emitters = None

        self.energy_observed = None
        self.spec_observed = None  # 1d array

    def _load_data(self):
        self.transit_data_dict = np.load(
            "data/data_local/transit/transit_data_dic_a%.3f_inc_%d.npz"
            % (self.spin, self.inc_angle),
            allow_pickle=True,
        )

    def SaveSpec(self, file_name):
        fname = "data/data_out/spectra_observed/" + file_name
        SaveXY(self.energy_observed, self.spec_observed, fname)

    def Transit(
        self,
        spectral_grid_emitters: SpectralDiskGrid,
        energy_observed,
    ):
        self.radi_disk_emitters = spectral_grid_emitters.radii
        self.energy_disk_emitters = spectral_grid_emitters.energies
        self.spec_disk_emitters = spectral_grid_emitters.spectra

        self.energy_observed = energy_observed
        self.spec_observed = self._ConvolSpec()

        return self.spec_observed

    def _ConvolSpec(self):

        radi_2d = self.transit_data_dict["radi_emitters"]
        garr_2d = self.transit_data_dict["redshift"]
        dSsc_2d = self.transit_data_dict["area_pixel_screen"]

        Neobs = len(self.energy_observed)
        spec_obs = np.zeros(Neobs)

        # Emin = np.min(EnSpec)
        # Emax = np.max(EnSpec)

        for i in range(Neobs):
            Ee = self.energy_observed[i] / garr_2d
            fSpec = reginter(
                (self.radi_disk_emitters, self.energy_disk_emitters),
                self.spec_disk_emitters,
                bounds_error=False,
                fill_value=0,
            )
            spec_new = fSpec((radi_2d, Ee))
            spec_arr = garr_2d**3 * spec_new * dSsc_2d
            spec_obs[i] = np.sum(spec_arr) / np.pi  # Because Inten_emis = Flux_emis/pi

        return spec_obs / self.distance**2
