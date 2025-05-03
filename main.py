import numpy as np
import sys
import traceback


from ziji_object import (
    BlackHole,
    AccretionDisk,
    Corona,
    Observer,
    ZijiParams,
)

from helper import (
    RunModel,
    InitLuminos,
    Error,
    Status,
    DBL_MIN
)

log_error = Error()
log_status = Status()



##Init Model
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



try:
    RunModel(ziji_params, folder_name)
    log_status.log_success(f"Success in {folder_name}")
except Exception as e:
    # Capture the full traceback as a string
    tb_str = traceback.format_exc()
    
    # Log the complete traceback and error message
    log_error.log(f"Error in {folder_name}: {e}\nTraceback: {tb_str}")
    log_status.log_failure(f"Failure in {folder_name}")