% Copyright (c) 2019, Christina Kaiser, Stefan Zeiske, Paul Meredith, Ardalan Armin, Swansea University
% All rights reserved.
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%__________________________________________________________________________
%   This program performs optical simulations for photovoltaic devices and was developed in MATLAB Version 9.4. 
%   It is aimed for determining the extinction coefficient
%   in the sub-gap spectral range, which is experimentally often not accessible. It requires as data input the device layer sequence, 
%   the layer thicknesses, the optical constants in the above-bandgap spectral range and the external
%   quantum efficiency spectra in the above and sub-gap spectral range.
%   For a detailed descripted and for citation please refer to: 
%   Determining Ultra-low Absorption Coefficients of Organic Semiconductors from the Sub-bandgap Photovoltaic External Quantum Efficiency.
%   Christina Kaiser, Stefan Zeiske, Paul Meredith, Ardalan Armin, Journal, xx, 2019.  https://doi.org/10.1002/adom.201901542
%__________________________________________________________________________
%   The function TMAT (and corresponding I_mat and L_mat) for calculating the optical field via the transfer matrix
%   method was taken from:
%   Accounting for Interference, Scattering, and Electrode Absorption to Make
%   Accurate Internal Quantum Efficiency Measurements in Organic and Other 
%   Thin Solar Cells, G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater., 22, 3293.
%_______________________________________________________________________
%   The numerical optimization problem was solved using the MATLAB(R) fminsearchbnd function as available
%   from https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon.
