import os, sys
import numpy as np
import multiprocessing
import pandas as pd
from scipy import constants
from astropy.io import ascii
import parse
from scipy import constants
from itertools import product
import tqdm

from IPython.core.debugger import Tracer

##############################################################################
# User inputs
nthreads = 3
drrecmode = "2"
zeta_vals = [f"{z:.2f}" for z in [0.05, 0.20, 0.40, 1.00, 2.00]]
logU_vals = [f"{u:.2f}" for u in np.arange(-4.0, -0.5, 0.5)]  # cm^-3
logPk_vals = [f"{p:.2f}" for p in np.arange(5.0, 10.0, 1.0)]

##############################################################################
# Delectronic recombination mode for MAPPINGS
with open ("/usr/local/share/mappings/data/switches.txt", "r") as f:
    data = f.readlines()
data[10] = f"{drrecmode}\n"
with open("/usr/local/share/mappings/data/switches.txt", "w") as f:
    f.writelines(data)

"""
Z/Z_sun     Default MAPPINGS abn file   linear_zeta_O_nebular   log10(Fe_free/Fe_tot)
0.05        GC_ZO_M1130                 0.086178462             -0.0737
0.20        GC_ZO_M0530                 0.331741775             -0.5347
0.40        GC_ZO_M0230                 0.605682622             -1.0593
1.00        GC_ZO_P0170                 1.342321468             -1.7548
2.00        GC_ZO_P0470                 2.53301829              -2.3095
"""

# Nebular abundance files
zeta_fnames = {
    "0.05": "GC_ZO_0086.abn",
    "0.20": "GC_ZO_0332.abn",
    "0.40": "GC_ZO_0606.abn",
    "1.00": "GC_ZO_1342.abn",
    "2.00": "GC_ZO_2533.abn"
}

dep_fnames = {
    "0.05": "Depletions_Fe-0.07.txt",
    "0.20": "Depletions_Fe-0.53.txt",
    "0.40": "Depletions_Fe-1.06.txt",
    "1.00": "Depletions_Fe-1.75.txt",
    "2.00": "Depletions_Fe-2.31.txt"
}

sb_fnames = {
    "0.05": "cont_a03t21isp_vm802.spectrum",
    "0.20": "cont_a03t22isp_vm802.spectrum",
    "0.40": "cont_a03t23isp_vm802.spectrum",
    "1.00": "cont_a03t24isp_vm802.spectrum",
    "2.00": "cont_a03t25isp_vm802.spectrum"
}

##############################################################################
# Wrapper function for calling Mappings
def map51(args):
    zeta_vals, logU, logPk = args

    # Create new subdirectory & cd into it
    os.system(f"[ -d logU{logU} ] || mkdir logU{logU}")
    os.chdir(f"logU{logU}")

    # Edit the .mv file
    with open(f"../phot-sph-tmp.mv", "r") as file:
        template_str = file.read()

    # Replace dummy values
    file_str = ""
    for ii, zeta in enumerate(zeta_vals):
        file_str += template_str.replace("NEBABUNDFNAME", zeta_fnames[zeta]).replace("NEBABUND", zeta).replace("DEPFNAME", dep_fnames[zeta]).replace("SBFNAME", sb_fnames[zeta])
        if ii < len(zeta_vals) - 1:
            file_str += "R\n"
        else:  
            file_str += "X\n"
    file_str = file_str.replace("LOGUVALUE", logU)
    file_str = file_str.replace("LPKVALUE", logPk)
    file_str = file_str.replace("DRRECMODE", drrecmode)

    # Save back to new file & delete template file
    mv_fname = f"phot-pp-logU{logU}-logPk{logPk}.mv"
    with open(mv_fname, "w") as file:
        file.write(file_str)

    # run mappings
    os.system(f"map51<{mv_fname} >> map_out.txt")

    # Move back into parent directory
    os.chdir("..")

    return

##############################################################################
# List of arguments to input to map51 wrapper
args_list = [[zeta_vals, logU, logPk] for logPk, logU in product(logPk_vals, logU_vals)]

##############################################################################
# For debugging
# map51(args_list[0])
# sys.exit(0)

# Before we run MAPPINGS, go through and remove old files
# for logU in logU_vals:
#     # Change into the directory
#     if os.path.exists(f"logU{logU}"):
#         os.chdir(f"logU{logU}")
#         os.system("rm *.csv")
#         os.system("rm *.ph6")
#         os.system("rm *.txt")
#         os.system("rm *.mv")
#         os.chdir("..")

# Run in multiple threads
print(f"Running HII region models...")
pool = multiprocessing.Pool(nthreads)
for _ in tqdm.tqdm(pool.imap_unordered(map51, args_list), total=len(args_list)):
    pass
