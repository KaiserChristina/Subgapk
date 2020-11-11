**Info**

This program performs optical simulations for photovoltaic devices and was developed in MATLAB Version 9.4. It is aimed for determining the extinction coefficient in the sub-gap spectral range, which is experimentally often not accessible. It requires as data input the device layer sequence, the layer thicknesses, the optical constants in the above-bandgap spectral range and the external quantum efficiency spectra in the above and sub-gap spectral range. For a detailed description and for citation please refer to: 
Kaiser, C.; Zeiske, S.; Meredith, P.; Armin, A. Determining Ultralow Absorption Coefficients of Organic Semiconductors from the Sub‐Bandgap Photovoltaic External Quantum Efficiency. Adv. Opt. Mater. 2019, 8 (1), 1901542. https://doi.org/10.1002/adom.201901542

The function TMAT (and corresponding I_mat and L_mat) for calculating the optical field via the transfer matrix
method was taken from:
Accounting for Interference, Scattering, and Electrode Absorption to Make
Accurate Internal Quantum Efficiency Measurements in Organic and Other 
Thin Solar Cells, G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater., 22, 3293.

The numerical optimization problem was solved using the MATLAB(R) fminsearchbnd function as available
from https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon.

__________________________________
**QuickStart**

To perform an exemplary simulation, download all files in one folder, open the file `<Subgapk.m>` and click run.
In the example, the photodiode consists of the layer sequence *glass/ITO/ZnO/material/MoO/Ag*. The optical constants
of each layer are provided in the file `<Optical constants.xls>` except for the sub-gap extinction coefficient of *material* that
is the photoactive layer of the photodiode. Moreover, three EQE spectra are provided corresponding to three photodiodes
with the above layer sequence, but different active layer thicknesses. 

To obtain the sub-gap extinction coefficient for a different material and layer sequence,
adapt the input section of the `<Subgapk.m>` file, add new optical constants to the `<Optical constants.xls>` 
and add new EQE spectra as .xls files to the folder. Change `<EQE_filename>` to your own files names in the input section of `<Subgapk.m>`.

Keep the format identical to the exemplary files, otherwise errors will be thrown.

