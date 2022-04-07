# bandas_comparadas
Python script to compare electronic band structure of the same system, or similar systems on the same k-point path. Example: comparing calculations with (without) spin-orbit interaction. Or for different functional/approximation results

Packages needed:

import os \
import numpy as np \
from mpl_toolkits.axes_grid1 import make_axes_locatable \
import matplotlib as mpl \
import matplotlib.pyplot as plt  
from scipy.interpolate import interp1d

Example for Sb2Te3, with SOC (red line) and without SOC (black dashed line).

![bandas_comparadas](https://user-images.githubusercontent.com/102771716/162301647-0f063e8d-3224-46a0-862b-e0d6518b9afc.png)

