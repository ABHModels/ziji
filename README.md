
# ziji – Simulating X-ray spectra of black hole binaries with returning radiation.

## Table of Contents
1. [Project Overview](#project-overview)
2. [Features](#features)
3. [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Installing](#installing)
4. [Usage](#usage)
5. [Folder Structure](#folder-structure)
6. [Contributing](#contributing)
7. [License](#license)
8. [Acknowledgements](#acknowledgements)

---

## Project Overview

`ziji` simulates the X-ray spectra of black hole X-ray binaries, incorporating the effects of returning radiation. It uses a disk+corona model, considering both thermal disk and reflection components. The model assumes a lamppost corona in Kerr spacetime and improves upon previous calculations by self-consistently determining the ionization parameter at each radial coordinate of the disk, resulting in more realistic reflection spectra.

References:    

[1] _Toward More Accurate Synthetic Reflection Spectra: Improving the Calculations of Returning Radiation_, Mirzaev et al., [Astrophys.J. 965 (2024](https://iopscience.iop.org/article/10.3847/1538-4357/ad303b)    
[2] _ZIJI: a model to calculate X-ray spectra of black hole X-ray binaries_, Mirzaev et al., [2406.01226](https://arxiv.org/abs/2406.01226)).


## Features

- Simulates X-ray spectra of black hole X-ray binaries with returning radiation effects.
- Incorporates a **disk+corona** model, considering both thermal disk and reflection components.
- Uses a **lamppost corona model** and assumes Kerr spacetime, with a **Novikov-Thorne disk**.
- The **reflection spectrum** is produced by both direct corona radiation and returning radiation of the thermal and reflection components, calculated with the actual spectrum illuminating the disk.
- Self-consistently calculates the **ionization parameter** at each radial coordinate of the accretion disk.


## Installation

### Prerequisites

- **Linux** system (compatible with all Linux platforms)
- **GCC 8+** (Supports C++17)
- **CMake 3.10+**
- **GSL** (GNU Scientific Library)
- **Python 3.8+**

### Installing

#### System Dependencies

- **For Ubuntu/Debian-based systems**:
   ```bash
   sudo apt update
   sudo apt install gcc g++ cmake libgsl-dev python3-pip
   ```

- **For CentOS/RedHat-based systems**:
   ```bash
   sudo yum install gcc gcc-c++ cmake gsl-devel python3-pip
   ```

#### Python Dependencies (using `requirements.txt`)

To install the required Python packages, use the provided `requirements.txt` file:

```bash
pip3 install -r requirements.txt
```

#### Universal Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/ABHModels/ziji.git
   cd ziji
   ```

2. **Compile the project**:
   ```bash
   cmake .
   make
   make cmake_clean
   chmod +x lib/reflector_linux
   ```


## Usage

### Running the Code

Running the project is simple. You can execute the main Python script as follows:
```bash
python3 main.py
```

### Changing Model Parameters
In the main.py script, users can modify various model parameters. For instance, in the following section of main.py:
```python
## Init Model
##################################

height = 5.
lambda_ = 0.1
rho_ = 100.

lumd, lumc = InitLuminos(rho_, lambda_)

ziji_params = ZijiParams(
    black_hole=BlackHole(mass=10, spin=0.998),
    accretion_disk=AccretionDisk(
        iron_abundance=1.0,
        scale_luminosity=lumd,
        color_factor=1.7,
        density_profile = np.full(100, 1e22),
    ),
    corona=Corona(gamma=2.0, ecut=300, height=height, scale_luminosity=lumc),
    observer=Observer(distance=1e6, inclination=45),
)

folder_name = f"folder_spectrum"
#######################################
```

Users can adjust parameters such as:

- `height`: The height of the corona.
- `lambda_`: Total luminosity of the source in Eddington units (λ = Ltot/LEdd).
- `rho_`: Relative luminosity of the corona (ρ = Lcorona/Ldisk).
- `black_hole`, `accretion_disk`, and `corona` properties, including mass, spin, iron abundance, and more.

All output data will be saved to the folder `ziji_disk/folder_spectrum`. Users can modify the name of this folder by changing `folder_spectrum` in the script.


## Folder Structure

Here is an overview of the important directories and files in the project:

```
/ziji
│
├── bin/                # Compiled binaries and shared libraries
├── data/               # Datasets for running simulations
├── external/           # External dependencies (e.g., xtensor, xtl)
├── ziji_disk/          # Output data from simulations
├── ziji_object/        # Python code for calculating integrations and finding spectra
├── zijiray/            # Main C++ codebase for ray tracing and model calculations
├── lib/                # Additional shared/static libraries, reflection models (e.g., xillver, reflionx, xstar)
├── README.md           # Project overview and installation guide
├── CMakeLists.txt      # CMake configuration file
├── requirements.txt    # Python dependencies
├── status_log.txt      # Log file for tracking statuses
├── error_log.txt       # Log file for tracking errors
├── main.py             # Main Python script for running the model
├── helper.py           # Helper Python script for utilities
├── data_manager.py     # Python script for managing data
└── ziji_param.json     # JSON file for storing model parameters
```

## Contributing

Contributions are welcome!
.......

## License

*(this section empty for now)*

## Acknowledgements

*(this section empty for now)*


