import os
from ctypes import *
import numpy as np

#######################
tmp = CDLL("bin/ray/libtransit_ray.so")

RayGen = tmp.RadiRedshArray

RayGen.argtypes = [c_double, c_double, c_wchar_p, c_int]
RayGen.restype = None

RscGen = tmp.rdiskLineScr

RscGen.argtypes = [c_int, c_double, c_double, c_wchar_p]
RscGen.restype = None
#########################


def dSscInit(Rsc, dph):
    Nradi, Nph = np.shape(Rsc)

    dSsc = np.zeros((Nradi, Nph))

    for i in range(Nph):
        r = Rsc[:, i]
        r2_, r1_ = np.zeros(Nradi), np.zeros(Nradi)
        r1_[0] = r[0]  # - (r[1] - r[0])/2 #change new vers 11 jan
        r2_[-1] = r[-1]  # + (r[-1] - r[-2])/2 #change new vers 11 jan

        for j in range(1, Nradi):
            r2_[j - 1] = (r[j] + r[j - 1]) / 2
            r1_[j] = r2_[j - 1]

        dSsc[:, i] = dph * (r2_**2 - r1_**2) / 2

    return dSsc


def InitArray(Nradi):
    phsc_ = np.load("data/data_cache/transit_ray/phsc.npy")
    rmax_ = np.load("data/data_cache/transit_ray/rscmax.npy")
    rmin_ = np.load("data/data_cache/transit_ray/rscmin.npy")

    Nph = len(phsc_)

    Psc, Rsc = np.zeros((Nradi, Nph)), np.zeros((Nradi, Nph))

    for i in range(Nph):
        Psc[:, i] = phsc_[i]
        LogArr = np.linspace(np.log10(rmin_[i]), np.log10(rmax_[i]), Nradi)
        Rsc[:, i] = 10**LogArr

    Xsc = Rsc * np.cos(Psc)
    Ysc = Rsc * np.sin(Psc)

    np.save("data/data_cache/transit_ray/Xsc.npy", Xsc)
    np.save("data/data_cache/transit_ray/Ysc.npy", Ysc)
    np.save("data/data_cache/transit_ray/Rsc.npy", Rsc)
    np.save("data/data_cache/transit_ray/Psc.npy", Psc)

    dSsc = dSscInit(Rsc, 2 * np.pi / Nph)
    np.save("data/data_cache/transit_ray/dSsc.npy", dSsc)


def check_file(file_path):

    if not os.path.exists(file_path):
        return False
    
    return True


def Save_data(spin, inc):
    dSsc = np.load("data/data_cache/transit_ray/dSsc.npy")
    # RadRed = np.load("data/RadiRedsh2d.npy")
    RadRed = np.load("data/data_cache/transit_ray/radi_redsh2d.npy")


    transit_data = {
        "radi_emitters": RadRed[:, :, 0],
        "redshift": RadRed[:, :, 1],
        "area_pixel_screen": dSsc,
    }


    fname = "data/data_local/transit/transit_data_dic_a%.3f_inc_%d.npz" % (spin, inc)
    np.savez(fname, **transit_data)




def TransitRay(Npar, spin, inc, cache = 0):
    fname = "data/data_local/transit/transit_data_dic_a%.3f_inc_%d.npz" % (spin, inc)
    check = check_file(fname)
    
    if check == False or cache==0:
        Nphsc = 200
        Nradi = 1000


        ph0 = -np.pi / 2 + 0.01
        ph1 = 2.0 * np.pi - np.pi / 2 + 0.01
        phsc = np.linspace(ph0, ph1, Nphsc, endpoint=False)

        # ph0 = -np.pi/2 + 0.01
        # ph1 = 2.*np.pi - np.pi/2 + 0.01 - 0.0001
        # phsc = np.linspace(ph0, ph1, Nphsc)
        np.save("data/data_cache/transit_ray/phsc.npy", phsc)

        RscGen(0, spin, inc, "data/data_cache/transit_ray/rscmax.npy")
        RscGen(1, spin, inc, "data/data_cache/transit_ray/rscmin.npy")

        InitArray(Nradi)

        # RayGen(spin, inc, "data/RadiRedsh2d.npy", Npar)
        RayGen(spin, inc, "data/data_cache/transit_ray/radi_redsh2d.npy", Npar)

        Save_data(spin, inc)



