# bandas_comparadas
Python script to compare electronic band structure of the same system, or similar systems on the same k-point path. Example: comparing calculations with (without) spin-orbit interaction. Or for different functional/approximation results

Packages needed:

import os \
import numpy as np \
from mpl_toolkits.axes_grid1 import make_axes_locatable \
import matplotlib as mpl \
import matplotlib.pyplot as plt \ 
from scipy.interpolate import interp1d \
