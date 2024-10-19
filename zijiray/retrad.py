import os
from ctypes import*
import numpy as np
import shutil


tmp = CDLL("bin/ray/libreturn_rad.so")


remred = tmp.RemRed

remred.argtypes = [c_double,c_double,c_wchar_p, c_int, c_int, c_int, POINTER(c_double)]  #(double spin, double incang, wchar_t* fname)
remred.restype = None


Nx = 800
Ny = 1200
# Nx = 100
# Ny = 200



def delete_and_recreate_folder(folder_path):
    # Check if the folder exists
    if os.path.exists(folder_path):
        # Remove the folder and its contents
        shutil.rmtree(folder_path)
    
    # Recreate the folder
    os.makedirs(folder_path)

def check_folder_and_file(folder_path, N, file_path):
    # Check if the folder exists

    if not os.path.exists(folder_path):
        return False
    
    # Count the number of files in the folder
    files_in_folder = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    if len(files_in_folder) != N:
        return False
    
    # Check if the specified file exists at the given location
    if not os.path.exists(file_path):
        return False
    
    return True


def DataGen(Npar, spin, robs_arr):

    delete_and_recreate_folder("data/data_local/returning_radiation/radrem_a%.3f"%spin)

    x0 =  np.arange(0, 6, dtype=float)
    length = len(x0)
    p_x0 = x0.ctypes.data_as(POINTER(c_double))

    paramarr = np.empty([0, length+1])

    for ind, rscr in enumerate(robs_arr):
        remred(spin,rscr,"data/data_local/returning_radiation/radrem_a%.3f/r%d.npy"%(spin,ind), Npar, Ny, Nx, p_x0)
        param = np.array(np.fromiter(p_x0, dtype=np.float64, count=length))
        param = np.append(rscr, param)
        paramarr = np.vstack([paramarr, param])

    with open("data/data_local/returning_radiation/robInPar_a%.3f.txt"%(spin), 'w') as f:
        # loop over each row in the array
        for row in paramarr:
            N = len(row)
            for i,val in enumerate(row):
                if(i != N-1):
                    f.write(str(val) + ' ')
                else:
                    f.write(str(val))
            f.write('\n')



def ReturningRay(Npar, spin, robs_arr, cache = 0):
    
    if cache == 0:
        DataGen(Npar, spin, robs_arr)
    else:
        folder_path = "data/data_local/returning_radiation/radrem_a%.3f"%spin
        N = len(robs_arr)
        file_path = "data/data_local/returning_radiation/robInPar_a%.3f.txt"%(spin)
        ch = check_folder_and_file(folder_path, N, file_path)

        if ch==False:
            DataGen(Npar, spin, robs_arr)

        