import json
from dataclasses import dataclass, asdict, field
from typing import Optional
import numpy as np


# Standalone filter function
def asdict_filter(items):
    return {k: v for k, v in items if not isinstance(v, np.ndarray)}


@dataclass
class BlackHole:
    spin: float  # Dimensionless
    mass: Optional[float] = None  # In Msun


@dataclass
class AccretionDisk:
    iron_abundance: float
    density: Optional[float] = None
    ionization: Optional[float] = None
    scale_luminosity: Optional[float] = None
    color_factor: Optional[float] = None

    density_profile: Optional[np.ndarray] = field(default=None, repr=False)
    ionization_profile: Optional[np.ndarray] = field(default=None, repr=False)


@dataclass
class Corona:
    gamma: float
    ecut: float
    height: float  # in rg
    scale_luminosity: Optional[float] = None


@dataclass
class Observer:
    distance: float  # in rg
    inclination: float  # in degrees


def format_large_numbers(obj):
    if isinstance(obj, dict):
        return {k: format_large_numbers(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [format_large_numbers(item) for item in obj]
    elif isinstance(obj, (int, float)):
        if abs(obj) >= 10000 or (0 < abs(obj) < 0.0001):
            return f"{obj:.1e}"
        else:
            return obj
    return obj


@dataclass
class ZijiParams:
    black_hole: BlackHole
    observer: Observer
    accretion_disk: Optional[AccretionDisk] = None
    corona: Optional[Corona] = None

    def save_to_json(self, file_path: str) -> None:
        # with open(file_path, "w") as file:
        #     raw_data = asdict(self, dict_factory=asdict_filter)
        #     formatted_data = format_large_numbers(raw_data)
        #     json.dump(formatted_data, file, indent=4)
        #     json_string = json.dumps(formatted_data, indent=4)
        #     file.write(json_string)
        raw_data = asdict(self, dict_factory=asdict_filter)
        json_string = self._to_json_string(raw_data, indent_level=0)
        with open(file_path, "w") as file:
            file.write(json_string)

    def _to_json_string(self, obj, indent_level):
        indent_space = "    "  # 4 spaces for indentation
        if isinstance(obj, dict):
            indent = indent_space * indent_level
            inner_indent = indent_space * (indent_level + 1)
            items = [
                f'\n{inner_indent}"{k}": {self._to_json_string(v, indent_level + 1)}'
                for k, v in obj.items()
            ]
            return f'{{{",".join(items)}\n{indent}}}'
        elif isinstance(obj, list):
            inner_indent = indent_space * (indent_level + 1)
            items = [
                f"\n{inner_indent}{self._to_json_string(v, indent_level + 1)}"
                for v in obj
            ]
            return f'[{",".join(items)}\n{indent_space * indent_level}]'
        elif isinstance(obj, float):
            # Apply custom scientific notation formatting
            if abs(obj) >= 1e4 or (0 < abs(obj) < 1e-3):
                return f"{obj:.1e}"
            else:
                return str(obj)
        elif obj is None:
            return "null"
        elif isinstance(obj, (int, str)):
            return json.dumps(obj)  # Use json.dumps for proper string and int handling
        return json.dumps(obj)

    @staticmethod
    def load_from_json(file_path: str) -> "ZijiParams":
        with open(file_path, "r") as file:
            data = json.load(file)

            black_hole_data = data.get("black_hole")
            accretion_disk_data = data.get("accretion_disk")
            corona_data = data.get("corona")
            observer_data = data.get("observer")

            black_hole = (
                BlackHole(**black_hole_data) if black_hole_data is not None else None
            )
            accretion_disk = (
                AccretionDisk(**accretion_disk_data)
                if accretion_disk_data is not None
                else None
            )
            corona = Corona(**corona_data) if corona_data is not None else None
            observer = Observer(**observer_data) if observer_data is not None else None

            return ZijiParams(
                black_hole=black_hole,
                accretion_disk=accretion_disk,
                corona=corona,
                observer=observer,
            )

    def __str__(self) -> str:
        return (
            f"ZijiParams:\n"
            f"  Black Hole: {self.black_hole}\n"
            f"  Accretion Disk: {self.accretion_disk}\n"
            f"  Corona: {self.corona}\n"
            f"  Observer: {self.observer}"
        )


class SpectralDiskGrid:
    def __init__(
        self, radii: np.ndarray, energies: np.ndarray, spectra: np.ndarray = None
    ):
        self.radii = np.copy(radii)
        self.energies = np.copy(energies)

        # If spectra is None, initialize it with a zeros array of the appropriate shape
        if spectra is None:
            self._spectra = np.zeros((len(self.radii), len(self.energies)))
        else:
            self._spectra = np.copy(spectra)  # Make a copy of the spectra if provided

    @property
    def spectra(self) -> np.ndarray:
        return self._spectra

    @spectra.setter
    def spectra(self, value: np.ndarray):
        if value.shape != (len(self.radii), len(self.energies)):
            raise ValueError(
                "Dimensions of spectra should match the length of radii and energies"
            )
        self._spectra = np.copy(value)  # Ensure a deep copy is made


if __name__ == "__main__":

    ziji_params = ZijiParams(
        black_hole=BlackHole(mass=10, spin=0.998),
        accretion_disk=AccretionDisk(
            iron_abundance=1.0,
            density=1e15,
            scale_luminosity=1e4,
            ionization=1000.0,
            ionization_profile=np.array([1, 2, 3]),
        ),
        corona=Corona(gamma=2.0, ecut=300, height=2.0, scale_luminosity=1e3),
        observer=Observer(distance=1e6, inclination=45),
    )

    fname_param = "ziji_param.json"
    ziji_params.save_to_json(fname_param)
    ziji_params_new = ZijiParams.load_from_json(fname_param)

    print(ziji_params_new)
