import os
from ctypes import *


tmp = CDLL("bin/ray/liblpg.so")

LPG = tmp.LPGeom

LPG.argtypes = [c_double, c_double]
LPG.restype = None


def LampostGeom(spin, hlp):
    LPG(spin, hlp)
