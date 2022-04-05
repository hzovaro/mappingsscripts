import os
import numpy as np
import pandas as pd
from astropy.io import ascii
import parse
from scipy import constants
from IPython.core.debugger import Tracer

##############################################################################
overwrite_hd5 = True
hd5_fname = f"HII_grid.hd5"
hd5_fname_rad = f"HII_grid_radial_profiles.hd5"

# Emission line fluxes to record
eline_lambdas_A = {
    "NEV3347": 3346.79,
    "OII3726": 3726.032,
    "OII3729": 3728.815,
    "NEIII3869": 3869.060,
    "HEI3889": 3889.0,
    "HEPSILON": 3970.072,
    "HDELTA": 4101.734,
    "HGAMMA": 4340.464,
    "HEI4471": 4471.479,
    "OIII4363": 4363.210,
    "HBETA" : 4861.325,
    "OIII4959": 4958.911,
    "OIII5007": 5006.843,
    "HEI5876": 5875.624,
    "OI6300": 6300.304,
    "SIII6312": 6312.060,
    "OI6364": 6363.776,
    "NII6548": 6548.04,
    "HALPHA": 6562.800,
    "NII6583": 6583.460,
    "SII6716": 6716.440,
    "SII6731": 6730.810,
    "SIII9069": 9068.600,
    "SIII9531": 9531.100
}
eline_list = eline_lambdas_A.keys()

##############################################################################
# Parse the output files
rows_list = []
rows_list_rad = []
logU_vals = [f.split("logU")[1] for f in os.listdir(".") if f.startswith("logU")]
for logU in logU_vals:
    # Change into the directory
    os.chdir(f"logU{logU}")

    # List of files of emission line intensities for each model
    spec_fnames = [f for f in os.listdir(".") if f.startswith("spec") and f.endswith(".csv")]
    rad_profile_fnames = [f for f in os.listdir(".") if f.startswith("phlss") and f.endswith(".ph6")]
    neb_structure_fnames = [f for f in os.listdir(".") if f.startswith("phsem") and f.endswith(".ph6")]
    run_numbers = [s.split("spec")[1].split(".csv")[0] for s in spec_fnames]

    # Fix dodgy columns
    os.system('find spec*.csv -type f -exec sed -i "" "s|, \\*| |g" {} \\;')

    # Store the emission line fluxes in a Pandas DataFrame
    for run_number in run_numbers:
        spec_fname = f"spec{run_number}.csv"
        rad_profile_fname = f"phlss{run_number}.ph6"
        neb_structure_fname = f"phsem{run_number}.ph6"
        atom_list = [fname.split(f"{run_number}")[0] for fname in [f for f in os.listdir(".") if f.endswith(f"{run_number}.csv") and not f.startswith("spec")]]
        row_dict = {}

        #######################################################################
        # Read the spec file to get the model parameters
        f = open(os.path.join(".", spec_fname))
        param_string = [f.readline() for ii in range(6)][-1].strip()
        f.close()

        # Read in the model parameters
        # print(param_string)
        print(f"In folder logU{logU}: {param_string}")
        try:
            zeta, logU_, logPk, drrecmode = parse.parse(
                'Run   :, HII,SPH,zeta={},logU={},logPk={},drrecmode={}',
                param_string)
        except:
            # This is only here in case we need to read old model files from before I put drrecmode in the model name.
            drrecmode = 2
            zeta, logU_, logPk = parse.parse(
                'Run   :, HII,SPH,zeta={},logU={},logPk={}',
                param_string)
        assert logU == logU_, "ERROR: logU in param_string != that in filename!"
        row_dict["Model"] = "HII"
        row_dict["Geometry"] = "SPH"
        row_dict["Nebular abundance (Z/Zsun)"] = float(zeta)
        row_dict["Abundance type"] = "GCZO"
        row_dict["log(U) (inner)"] = float(logU)
        row_dict["log(Q) (inner)"] = float(logU) + np.log10(constants.c * 1e2)
        row_dict["log(P/k)"] = float(logPk)
        row_dict["Stellar abundance (Z/Zsun)"] = float(zeta)
        row_dict["IMF"] = "Salpeter"
        row_dict["drrecmode"] = drrecmode

        #######################################################################
        # Read the spec file to get the line intensities
        num_lines = sum(1 for line in open(os.path.join(".", spec_fname)))
        table = ascii.read(os.path.join(".", spec_fname),
            data_start=-(num_lines - 56), header_start=-(num_lines - 54))
        for eline in eline_list:
            eline_idx = np.argmin(
                np.abs(table['Lambda(A)'] - eline_lambdas_A[eline]))
            row_dict[eline] = table['Flux (HB=1.0)'][eline_idx]

        #######################################################################
        # Get the model's radius and temperature values
        # NOTE: we get R, T from this file rather than the ion files because we get better precision this way.
        # ALSO NOTE: this doesn't parse most of the columns correctly!! But R and T are OK.
        d = pd.read_csv(rad_profile_fname, skiprows=34, delim_whitespace=True)
        row_dict_rad = row_dict.copy()
        row_dict_rad["R (cm)"] = d["Dist."].values
        row_dict_rad["T_e (K)"] = d["Te"].values

        #######################################################################
        # Read the phsem file to get the nebular structure
        d = pd.read_csv(neb_structure_fname, skiprows=37, skipfooter=3)
        # Fix the bad formatting
        l = [s.strip().split(" ")[1] for s in d.iloc[:, 0]]
        l[0] = 0.0
        d.iloc[:, 0] = [float(s) for s in l]
        colnames = d.columns[1:]
        for cc in range(len(d.columns)):
            d = d.rename(columns={d.columns[cc]: f"col {cc}"})
        for cc in range(len(colnames)):
            d = d.rename(columns={f"col {cc}": colnames[cc].strip()})
        # Trim remaining 2 columns
        d = d.iloc[:, :-2]
        row_dict_rad["n_H (cm^-3)"] = d["[6]nH"].values
        row_dict_rad["n_e (cm^-3)"] = d["[7]ne"].values
        row_dict_rad["n_tot (cm^-3)"] = d["[8]nt"].values
        row_dict_rad["SII6716 (radial values)"] = d["[19]6716"].values
        row_dict_rad["SII6731 (radial values)"] = d["[18]6731"].values
        row_dict_rad["HBETA (radial values)"] = d["[11]HBeta"].values

        #######################################################################
        # Read ions, if the files exist
        if len(atom_list) > 0:
            # For each element, read the ionisation state file
            for atom in atom_list:
                d = pd.read_csv(f"{atom}{run_number}.csv", skiprows=2)
                for col in d.columns:
                    d = d.rename(columns={col: col.strip()})

                # Iterate through each ionisation state, add to row
                for c in d.columns[11:-1]:
                    ion_state = c.split(" ")[1]
                    row_dict_rad[f"{atom} {ion_state}"] = d[c].values

            # Append to list
            rows_list.append(row_dict)
            rows_list_rad.append(row_dict_rad)

    #######################################################################
    # Move back into parent directory
    os.chdir("..")

# Convert to DataFrame
df = pd.DataFrame(rows_list)
df_rad = pd.DataFrame(rows_list_rad)

# Save to hd5 file
if overwrite_hd5:
    df.to_hdf(hd5_fname, key=f"HII")
    df_rad.to_hdf(hd5_fname_rad, key=f"HII")
else:
    # If HDF file does not exist, then create it
    if os.path.exists(hd5_fname):
        df_old = pd.read_hdf(hd5_fname)
        df = df_old.append(df, ignore_index=True)
    df.to_hdf(hd5_fname, key=f"HII")

    if os.path.exists(hd5_fname_rad):
        df_old = pd.read_hdf(hd5_fname_rad)
        df_rad = df_old.append(df_rad, ignore_index=True)
    df_rad.to_hdf(hd5_fname_rad, key=f"HII")

