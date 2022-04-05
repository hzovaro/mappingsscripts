from __future__ import division, print_function
from collections import OrderedDict as OD
from datetime import datetime  # For current date

import numpy as np
import pandas as pd


"""
Generate dust depletion input files for MAPPINGS.

The depletions are parametrised by DFe = log10(Fe_gas/Fe_total).

We optionally parametrise DFe by the oxygen abundance (very roughly).  So
given an input oxygen abundance (12 + log O/H), we can write out a MAPPINGS
input file with hopefully appropriate depletions.

Adam D. Thomas 2020 May
"""
__version__ = 0.1




##########################################
########### Inputs begin here ############
##########################################

output_dir = ("/Users/adamt/Documents/Work_local/Data/MAPPINGS_grids/"
              "MAPPINGS_grid_11_HII/custom_inputs/depletions/")
output_filename_template = "Depletions_Fe{0:.2f}.txt"
# List the DFe = log10(Fe_gas/Fe_total) values to write outputs for:
list_DFe = [0.00]  # -2.25, -1.00, 0.00, -1.5
# List 12 + log O/H values to use to estimate DFe and also write outputs for:
list_12OH = [7.297, 7.695, 7.996, 8.281, 8.434, 8.542, 8.695, 8.804, 8.888,
             8.957, 9.015, 9.093, 9.164, 9.261, 9.340, 9.407]
# Set to True to also write all outputs as csv files (a more friendly format):
write_csvs_too = False

##########################################
############ Inputs end here #############
##########################################




# ADT: I've reverse-engineered this table of Mike Dopita's interpolation
# scheme from the depletion files in mappingsv/lab/abund/unified_depletion/.
# The gradient T and intercept S data are for the scaling relation:
#    log10(X_gas/X_tot) = min(0, T * DFe + S)
colnames = ["Element_number", "Element_name", "Gradient_T", "Intercept_S"]
DFe_scaling_table = [
                        (1 ,       "H" ,       0       ,     0.0     ),
                        (2 ,       "He",       0       ,     0.0     ),
                        (3 ,       "Li",       0.47    ,     0.49    ),
                        (4 ,       "Be",       0.466667,     0.303333),
                        (5 ,       "B" ,       0.468571,     0.127143),
                        (6 ,       "C" ,       0.075556,    -0.038889),
                        (7 ,       "N" ,       0.08    ,     0.08    ),
                        (8 ,       "O" ,       0.176   ,     0.158667),
                        (9 ,       "F" ,       0.96    ,     1.35    ),
                        (10,       "Ne",       0       ,     0.0     ),
                        (11,       "Na",       0.466667,     0.283333),
                        (12,       "Mg",       0.78    ,     0.476667),
                        (13,       "Al",       0.78    ,     0.476667),
                        (14,       "Si",       0.886667,     0.623333),
                        (15,       "P" ,       0.733333,     0.99    ),
                        (16,       "S" ,       0       ,     0.0     ),
                        (17,       "Cl",       0.96    ,     1.35    ),
                        (18,       "Ar",       0       ,     0.0     ),
                        (19,       "K" ,       0.47    ,     0.09    ),
                        (20,       "Ca",       1.594286,     0.438571),
                        (21,       "Sc",       0.78    ,     0.483333),
                        (22,       "Ti",       1.594286,     0.438571),
                        (23,       "V" ,       1.948571,     0.757143),
                        (24,       "Cr",       1.125   ,     0.24    ),
                        (25,       "Mn",       0.666667,    -0.273333),
                        (26,       "Fe",       1       ,     0.0     ),
                        (27,       "Co",       1.05    ,    -0.06    ),
                        (28,       "Ni",       1.16    ,     0.17    ),
                        (29,       "Cu",       0.551111,    -0.074444),
                        (30,       "Zn",       0.48    ,     0.52    ),
]
table_as_cols = list(zip(*DFe_scaling_table))
DF_DFe = pd.DataFrame(OD(zip(colnames, table_as_cols)))
DF_DFe.set_index("Element_number", inplace=True)
# Sanity checks: we must include all elements from 1 to 30, depletions must
# increase with DFe, there is trivial scaling for Fe itself.
assert DF_DFe.index.tolist() == list(range(1,31))
assert np.all(DF_DFe["Gradient_T"].values >= 0)
assert DF_DFe.at[26,"Gradient_T"] == 1.0
assert DF_DFe.at[26,"Intercept_S"] == 0.0


file_header_template = """\
%
%  FORMAT
%  Title
%  # entries [type] (I2,x,I2)
%  Z  Offset/Scale
%  .
%  .
%  .
%
%  Note:
%  type = 0, or no type present (ie old files)
%        negative/0.0 offsets/scales are logarithmic, ie 0.00 is no change
%        positive offsets/scales are linear, ie 1.00 is no change
%  type = 1:
%        All offsets/scales are logarithmic, ie +0.00 is no change
%  type = 2:
%        All offsets/scales are linear , ie x1.00 is no change
%  type = 3:
%        All offsets/scales are log+12 , ie H = 12.00
%
%  Elements not included in the data file map.prefs will not be read
%
% Depletion file for Jenkins (2014) data with MAD interpolation
% Written by python script v{0} on {1}
"""


def calculate_depletions(DFe):
    """
    Calculate the depletion of each element.

    Where DFe = log10(Fe_free/Fe_tot), calculates the depletions for all
    other elements, using the formula:
        log10(X_gas/X_tot) = min(0, T * DFe + S)
    where the gradient T and intercept S are from the hard-coded table above.

    Outputs a table with a row for each element and a column for depletions.
    """
    if DFe > 0:
        raise ValueError("DFe can't be positive: DFe = log10(Fe_gas / Fe_tot)")

    DF_out = DF_DFe.copy()
    T = DF_out["Gradient_T"].values
    S = DF_out["Intercept_S"].values
    DF_out["log10(X_gas/X_tot)"] = np.minimum(0.0, T * DFe + S)
    return DF_out[["Element_name", "log10(X_gas/X_tot)"]]



def write_depletion_file(DF_out, DFe, extra_text=""):
    """
    Write out a depletion file for MAPPINGS for a particular DFe value.
    DF_out is the table of depletions, as a pandas DataFrame.
    extra_text is appended to the file header.
    """
    print("Saving depletion file for DFe = {0:.4f} ...".format(DFe))
    # Sanity check
    assert DF_out.at[26,"log10(X_gas/X_tot)"] == DFe
    assert isinstance(extra_text, str)

    fname = output_dir + output_filename_template.format(DFe)
    with open(fname, "w") as f:
        today = datetime.date(datetime.now())
        header1 = file_header_template.format(__version__, today)
        f.write(header1)
        if extra_text:
            f.write("% " + extra_text + "\n")
        title = "Base Depletion: log10 Fe_gas/Fe_tot = {0:0.4f}\n30  1"
        f.write(title.format(DFe) + "\n")
        elems, deps = DF_out.index.values, DF_out["log10(X_gas/X_tot)"].values
        to_add = []
        for elem, dep in zip(elems, deps):
            to_add.append("{0:>2d}  {1:>5.2f}".format(elem, dep))
        f.write("\n".join(to_add) + "\n")

    if write_csvs_too:
        # Also save as a simple csv for convenience
        print("Saving csv depletion file for DFe =", DFe, "...")
        DF_out.to_csv(fname[:-4] + ".csv")



# Table of coefficients for approximate conversion of 12 + log O/H -> DFe
line_table_colnames = ["Line_number",  "Gradient", "Intercept"]
DFe_scaling_table = [
                           (1 ,         -0.07,       0.4751 ),
                           (2 ,         -0.43,       3.2351 ),
                           (3 ,         -2.01,      16.1101 ),]
line_table_as_cols = list(zip(*DFe_scaling_table))
DF_convert = pd.DataFrame(OD(zip(line_table_colnames, line_table_as_cols)))
DF_convert.set_index("Line_number", inplace=True)
assert DF_convert.index.tolist() == list(range(1,4))
assert np.all(DF_convert["Gradient"].values < 0)
assert np.all(DF_convert["Intercept"].values > 0)



def convert_O_abundance_to_DFe(log_OH_12):
    """
    Use an input oxygen abundance to estimate a representative Fe depletion.
    (There are issues with even trying to do this, but I'm attempting to).
    Use a 3-part piecewise linear function designed in particular to give:
        12 + log O/H     OH/OH_Sun     Z/ZSun      DFe
            8.157          0.25         0.205     -0.25
            8.511          0.56         0.501     -1.00
            8.760          1.00         1.000     -1.50
    The function changes slope twice, at and below the lowest of these points.
    The required coefficients are hard-coded above.

    input:  12 + log O/H value as a float
    output: Estimate for DFe = log10(Fe_gas/Fe_tot)
    """
    if (log_OH_12 < 6.0) or (log_OH_12 > 10.0):
        raise ValueError("Bad 12 + log O/H value: {0:.4f}".format(log_OH_12))

    # y = mx+b
    x = log_OH_12
    m1, m2, m3 = DF_convert["Gradient"].tolist()
    b1, b2, b3 = DF_convert["Intercept"].tolist()
    DFe = min(0.0, m1*x + b1, m2*x + b2, m3*x + b3)
    return DFe



# Iterate over user inputs, writing requested files
for DFe_i in list_DFe:
    DF_i = calculate_depletions(DFe_i)
    write_depletion_file(DF_i, DFe_i)

for log_OH_12 in list_12OH:
    DFe_i = convert_O_abundance_to_DFe(log_OH_12)
    print("Note: 12 + log O/H = {0:.4f}".format(log_OH_12),
          "converted to DFe = {0:.4f}".format(DFe_i))
    DF_i = calculate_depletions(DFe_i)
    note = "Base Fe depletion guesstimated from 12 + log O/H = {0:.4f}"
    write_depletion_file(DF_i, DFe_i, extra_text=note.format(log_OH_12))


print("Done.")
