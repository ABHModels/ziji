import sys
import logging
import numpy as np
import os

from ziji_object import (
    ZijiParams,
    Iterator,
)
from ziji_object import GenGlobalRadiArray
from data_manager import manage_files
from zijiray import TransitRay
from zijiray import LampostGeom
from zijiray import ReturningRay

def GenRelativity(ziji_params: ZijiParams, Npar = 100):
    spin = ziji_params.black_hole.spin
    inc = ziji_params.observer.inclination
    hlp = ziji_params.corona.height
    TransitRay(Npar=Npar, spin=spin, inc=inc, cache=1)
    LampostGeom(spin=spin, hlp=hlp)
    robs_arr = GenGlobalRadiArray(spin, N=100)
    ReturningRay(Npar=Npar, spin=spin, robs_arr = robs_arr, cache = 1)


def RunModel(ziji_params, folder_name):
    GenRelativity(ziji_params)
    iterator = Iterator(ziji_params, mode=1)
    Xi_ion = iterator.InitModel()
    iterator.Loop()
    ziji_params.accretion_disk.ionization = Xi_ion

    fname_param = "ziji_param.json"
    ziji_params.save_to_json(fname_param)

    source_directory_path = "data/data_out/"
    target_base_directory = "ziji_disk"
    new_folder_name = folder_name
    manage_files(source_directory_path, target_base_directory, new_folder_name)

def InitLuminos(rho_, lambda_):
    lumbd = lambda_ / (1.0 + rho_)
    lumcor = rho_ * lumbd

    return lumbd, lumcor


def DeleteAll():
    os.system('find data -type f -not -path "data/data_local/global_data/*" -delete')



class Error:
    def __init__(self, file_error="error_log.txt"):
        # Create a custom logger for Error
        self.logger = logging.getLogger('ErrorLogger')
        self.logger.setLevel(logging.ERROR)
        
        # Prevent adding multiple handlers if logger already has them
        if not self.logger.handlers:
            # Create handlers
            fh = logging.FileHandler(file_error)
            fh.setLevel(logging.ERROR)
            
            # Create formatter and add it to handlers
            formatter = logging.Formatter('%(asctime)s - %(message)s')
            fh.setFormatter(formatter)
            
            # Add handlers to the logger
            self.logger.addHandler(fh)

    def log(self, txt):
        self.logger.error(txt)


class Status:
    def __init__(self, file_status="status_log.txt"):
        # Create a custom logger for Status
        self.logger = logging.getLogger('StatusLogger')
        self.logger.setLevel(logging.INFO)
        
        # Prevent adding multiple handlers if logger already has them
        if not self.logger.handlers:
            # Create handlers
            fh = logging.FileHandler(file_status)
            fh.setLevel(logging.INFO)
            
            # Create formatter and add it to handlers
            formatter = logging.Formatter('%(asctime)s - %(message)s')
            fh.setFormatter(formatter)
            
            # Add handlers to the logger
            self.logger.addHandler(fh)

    def log_success(self, txt="Program finished successfully"):
        self.logger.info(txt)

    def log_failure(self, txt="Program did not finish successfully"):
        self.logger.info(txt)




DBL_MIN = sys.float_info.min