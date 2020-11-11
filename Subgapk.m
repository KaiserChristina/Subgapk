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
 
function Numerical_k_from_EQE
close all
clc
clear
warning('off')
format shortEng
 
%=============define global variables======================
global layers
global lambda_Vis
global activeLayer
global thicknesses
global EQE
global n 
global wavelength
global t_active
global A_e_final
global alpha
global t_final
 
%=================user input=========================
% The optical constants of the device layers must be listed as shown in the
% Optical constants.xls file. The desired wavelength range (lambda_start to
% lambda_end) should not exceed the wavelength span for any material given in this file, 
% except for the sub-bandgap extinction coefficient of the active layer
% that is yet to be determined. The real part of the complex refractive
% index of the active layer must span the wavelength range
% lambda_start to lambda_end, whereas the imaginary part of the complex
% refractive index (extinction coefficient) must span lambda_start to lambda_bandgap. 
 
nk_filename = 'Optical constants.xls'; 
                
EQE_filename = {...                     % EQE spectra need to have the same wavelength range.
    'EQE_50nm.xls'
    'EQE_100nm.xls'
    'EQE_150nm.xls'
   };
 
lambda_step = 5;                        % stepsize within desired wavelength range
lambda_start = 400;                     % minimum of wavelength range
lambda_bandgap = 600;                   % maximum of the wavelength range in which the experimental extinction coefficient of the active layer is known
lambda_end = 1200;                      % maximum of wavelength range for all other experimental optical constants
 
layers = {'glass' 'ITO' 'ZnO' ...
    'material' 'MoO' 'Ag'};             % material layer sequence in the device starting from light incident side
thicknesses = [0 104 37 NaN 7 130];     % value of NaN is specified by t_active
activeLayer = 4;                        % index of material layer where photocurrent is generated
t_active = [53 108 150];                % experimental activeLayer thickness in nm; requires as many entries as length of EQE_filename
 
t_iterations = 20;                      % number of iterations for numerically obtaining t
k_iterations = 30;                      % number of iterations for numerically obtaining k
t_interval = '10';                      % interval (+/- x nm) in which the active layer thickness is optimized; set to 0 if t_active is precise
k_interval = 'k/10';                    % interval around previous k (0 to k + k_interval) in which k is optimized 

% The numerical optimization of k below the bandgap energy CAN have converge to extremely high
% values (k > 100) if k decreases sharply at the bandgap edge.
% Therefore, the numerical k is not allowed to exceed k (at the bandgap) for 
% wavelengths above the bandgap and below bandgap+stepsize*5.
% If lambda_bandgap does not correspond to the bandgap energy (or the decrease in k at the bandgap is slow), it is
% advised to select k_fixed = 'no'. It that case, the upper constrained for the numerical k above and below the
% bandgap are the same.

k_fixed = 'yes'; 

%=================End of User input=========================
 
lambda_Vis = lambda_start:lambda_step:lambda_bandgap;       % wavelength range where experimental k of the active layer is available; IQE is determined here
wavelength = lambda_start:lambda_step:lambda_end;           % full wavelength range over which k can be fit based on IQE and EQE
 
%----Load optical constants and EQE, and check for errors in the user input
if length(EQE_filename) ~= length(t_active)
    disp('Error occured: Number of EQE_filename does not match the size of t_active.')
    return
elseif lambda_bandgap > lambda_end 
    disp(['Error occured: Wavelength range of the experimental extinction coefficient' newline ...
        'of the active layer is larger than the wavelength range over which the extinction coefficient should be fit.'])
    return
elseif length(thicknesses) ~= length(layers)
    disp('Error occured: Each layer in layers requires a thickness.')
    return
end
 
try                                                         % executes a piece of code
    n = zeros(length(layers),length(wavelength));
    for index = 1:length(layers)
        n(index,:) = LoadRefrIndex(layers{index},wavelength, nk_filename); % complex refractive index n is interpolated depending on lambda
    end
catch                                                       % shows warning if data readout from excel is file not possible
    disp(['Error occured during readout of the optical constants.' newline 'Layers must correspond to the name in the "Optical_constants.xls" file.'])
    return
end
 
try
    EQE = LoadFile(EQE_filename,wavelength);
catch
    disp(['Error occured during EQE readout. Check for correct EQE spectra filename and location.' newline 'Wavelength range specified by lambda_start and lambda_end must be within the experimental EQE spectral range.' newline ...
        'Maximum and minimum wavelength of all EQE spectra must be the same.'])
    return
end
 
%-----------------Thickness pre-optimization-------------------------------
                
opts = optimset('fminsearch');                      % optimization option structure
opts.Display = 'off';
opts.TolX = 1e-3;                                   % Limits precision of active layer determination
opts.MaxFunEvals = t_iterations; 
 
A_s = zeros(length(lambda_Vis),length(t_active)); 
t_final = zeros(length(t_active),1); 
for index = 1:length(t_active)                      % index of device of which the thickness is optimized
    global EQE_device                               % define global variable EQE_device for the first objective function thickness_opt that optimizes the ative layer thickness of ONE device
    EQE_device = EQE(1:length(lambda_Vis),index);   % restrict the thickness optimization to the wavelength range, where optical constants of all layers are known
    x0 = [t_active(index)];                         % starting point for the thickness optimization is the experimental thickness
    
    % thickness optimization using the Nelder-Mead simplex algorithm; x=FMINSEARCHBND(fun,x0,LB,UB,options)
    % thickness is allowed to vary within +/- t_interval
    x = fminsearchbnd(@thickness_opt,x0,[t_active(index)- eval(t_interval)],[t_active(index)+eval(t_interval)], opts); 
    t_final(index,:) = round(x(1));  
    
    method = 1;                                     % Method 1 performs a standard transfer matrix simulation over a broad wavelength range.
    
    A_s(:,index) = TMAT(lambda_Vis,[ ], n(:,1:length(lambda_Vis)), thicknesses, t_final(index) , method , activeLayer); % The active layer absorption (A_s) is simulated by the transfer matrix method using the optimized thickness.

    IQE = EQE(1:length(lambda_Vis),:)./ A_s;        % The IQE is estimated in the visible spectral range, where the optical constants of all materials are known.
    
    % Under the assumption that the IQE is excitation energy independent, the absorption (A = EQE*IQE_(mean)^(-1)) is
    % estimated for all devices over the spectral range of the experimental EQE spectra.
    % A_e_final is a matrix containing the estimated absorption of all devices over the full spectral range; size(A_e_final) = length(wavelength) x length(t_active)
    A_e_final = bsxfun(@rdivide, EQE(1:length(wavelength),:), mean(IQE));
end
    t_final
%---------------Plotting the active layer absorption
 
% This figure shows the simulated active layer absorption in the visible
% spectral using the transfer matrix method.
 
figure('Name',strcat(layers{activeLayer}," ",'absorption from transfer matrix simulation','_',datestr(now, 'ddmmyy_HHMM')),'NumberTitle','off','DefaultAxesFontSize',8)
plotString = '';
for index = 1:length(t_active)
    plotString = strcat(plotString, ['lambda_Vis , A_s(:,', num2str(index), '),']);
end
eval(['plot(',plotString,'[],[],''LineWidth'',2)'])
xlabel('Wavelength (nm)');
ylabel('Absorption');
set(gcf,'units','centimeters','position',[0,0,11,8])
legend(strcat(string(t_final),' nm'), 'Location', 'best', 'Interpreter','none');
drawnow
 
% This figure depicts the IQE of the device in the visible spectral range.
% In the case of significant wavelength dependence (IQE deviates
% from the mean IQE by more than 30 %), it is not advised to use this
% method for obtaining sub-bandgap absorption. The deviation can be
% minimized by improving the accuracy of the model (e.g. precise optical
% constants and layer thicknesses) or by increasing t_interval within reason to allow for
% an optimization of the active layer thickness over a wider range. However, the deviation 
% depends mostly on the accuracy of the active layer optical constants in the visible
% spectral range.
 
figure('Name',strcat('Device IQE','_',datestr(now, 'ddmmyy_HHMM')),'NumberTitle','off','DefaultAxesFontSize',8)
plotString = '';
for index = 1:length(t_active)
    plotString = strcat(plotString, ['lambda_Vis , IQE(:,', num2str(index), '),']);
end
eval(['plot(',plotString,'[],[],''LineWidth'',2)'])
xlabel('Wavelength (nm)');
ylim([-0.05 1.05])
ylabel('IQE');
set(gcf,'units','centimeters','position',[0,0,11,8])
legend(strcat(string(t_final),' nm'), 'Location', 'best', 'Interpreter','none');
drawnow

 
% This figure shows the estimated active layer absorption (A_e = EQE/IQE_mean)
% over the full wavelength range assuming that the IQE is excitation energy independent.
% A_e and A_s are the same, if the IQE is 1.
 
figure('Name',strcat(layers{activeLayer}," ",'absorption','_',datestr(now, 'ddmmyy_HHMM')),'NumberTitle','off','DefaultAxesFontSize',8)
plotString = '';
for index = 1:length(t_active)
    plotString = strcat(plotString, ['wavelength , A_e_final(:,', num2str(index), '),']);
end
eval(['plot(',plotString,'[],[],''LineWidth'',2)'])
xlabel('Wavelength (nm)');
ylabel('Absorption');
set(gcf,'units','centimeters','position',[0,0,11,8])
legend(strcat(string(t_final),' nm'), 'Location', 'best', 'Interpreter','none');
drawnow
 
%---------------k optimization
opts = optimset('fminsearch');              %  optimization option structure
opts.Display = 'off';
opts.TolX = 1e-8; % precision here must be lower than the minimum k
opts.MaxFunEvals = k_iterations;
 
k_numerical = zeros(length(wavelength),1);
figure('Name',strcat('Numerical k of'," ", layers{activeLayer},'_',datestr(now, 'ddmmyy_HHMM')),'NumberTitle','off','DefaultAxesFontSize',8)
for index = 1:length(wavelength)
    alpha = wavelength(index); % defines wavelength at which optimization is run within the second objective function k_opt
%     
%   fprintf('Current Wavelength %d \n', alpha); % updates status of the simulation
%     
    % The numerical optimization does not require the experimental
    % extinction coefficient once the IQE in the visible spectral range is
    % determined. However, the experimental k at the lowest wavelength is chosen 
    % as a starting point for the simulation of k at higher wavelengths.
  
    if index == 1                   % the initial k is approximated from first experimental k
        k = imag(n(activeLayer,1)); 
    end
    x0 = k;                         % starting point of the optimization
    
    % The numerical optimization of k below the bandgap energy CAN have converge to extremely high
    % values (k > 100) if k decreases sharply at the bandgap edge.
    % Therefore, the numerical k is not allowed to exceed k at the bandgap for wavelengths above the bandgap and below bandgap+stepsize*5.
    % Change the value "5" to increase/decrease the wavelength range in which k is fixed. 
    
    if strcmp(k_fixed,'yes') && index > find(wavelength == lambda_bandgap) && index < find(wavelength == lambda_bandgap) + 5
    x = fminsearchbnd(@k_opt,x0,[],[k], opts);  % x = FMINSEARCHBND(fun,x0,LB,UB,options)
    else
    x = fminsearchbnd(@k_opt,x0,[],[k + eval(k_interval)], opts);   
    end
    k = x(1);                       % defining k for the next loop
    k_numerical(index) = k;         % compiling k
    
    %----plotting k------------------
    semilogy(wavelength(1):lambda_step:alpha, k_numerical(1:index),'r-','LineWidth',2) % plotting the numerical k in real-time
    hold on
    xlabel('Wavelength (nm)');
    ylabel('Extinction coefficient');
    set(gcf,'units','centimeters','position',[5,5,11,8])
    drawnow 
end
    assignin('base','k', k_numerical) % defining variable for export
    
 
function f = thickness_opt(x) % pre-optimizing the active layer thickness
global activeLayer
global n
global lambda_Vis
global thicknesses
global EQE_device
% redefining thicknesses to loose the global variable property in order to substitute for the active layer thickness depending on the device
t = thicknesses; 
activeLayer_thickness = x(1); % defining the variable
t(activeLayer) = activeLayer_thickness;
method = 1;    
% Method 1 executes the standard transfer matrix method resulting in the
% active layer absorption from the experimental optical constants.
% k = [] because k will be loaded from experimental data
A_s = TMAT(lambda_Vis,[ ], n(:,1:length(lambda_Vis)), t, activeLayer_thickness, method, activeLayer);% transfer matrix method
IQE = EQE_device(1:length(lambda_Vis),1)./A_s.'; 
A_e = EQE_device(1:length(lambda_Vis),1)./mean(IQE);
% The objective function f minimizes the relative difference between the
% simulated and the estimated absorption. f has to be a scalar and
% therefore the difference is summed up over all wavelengths. The absolute
% difference is normalized to the estimated absorption in order to minimize
% difference at all wavelengths independent of the absolute absorption value. 
f = sum(((A_s.' - A_e).^2)./A_e); %
 
function f = k_opt(x)       % numerical optimization of k
global activeLayer
global t_active
global n
global wavelength
global thicknesses
global alpha                % one wavelength
global t_final              % optimized thickness
global A_e_final 
 
A_e = A_e_final((find(wavelength == alpha)),:).'; % estimated absorption (A = EQE*IQE_(mean)^(-1)) at one wavelength for all devices
k = x(1);                   % defining the variable
method = 2;                 % Method 2 executes the transfer matrix method with a fixed
                            % thickness (t_final) and at ONE wavelength. 

A_n = zeros(length(wavelength),length(t_active)); 
for index = 1:length(t_active) % indexes over all devices
    A_n(:, index) = TMAT(wavelength, k , n , thicknesses, t_final(index) , method , activeLayer); % transfer matrix method
end
% find numerical absorption of all devices corresponding to the wavelength at which k is currently optimized
A_n = A_n((find(wavelength == alpha)),:).';
% The objective function f minimizes the relative difference between the
% estimated and the numerical absorption. f must be a scalar and
% therefore, the difference is summed up over all wavelengths. The absolute
% difference is normalized to the estimated absorption in order to minimize
% difference at all wavelengths independent of the absolute absorption value. 
f = sum(((A_e - A_n).^2)./A_e); % normalized objective function

 
function Absorption_activeLayer = TMAT(wellenlange, k, n, thicknesses, activeLayer_thickness, method, activeLayer)
 
thicknesses(activeLayer) = activeLayer_thickness; 
t = thicknesses;
 
%--------------Method 1 & 2------------------------------------

if method == 1                                  % Method 1: Standard transfer matrix simulation over a wavelength range. k is not varied.
    n = n(:,1:length(wellenlange));                                     
elseif method == 2                              % Method 2: Transfer matrix simulation for one wavelength only. 
                                                % The extinction coefficient of the active layer is substituted at one wavelength. 
                                                
    global alpha                                % call the wavelength at which k is currently optimized 
                                        
    n(activeLayer, find(wellenlange == alpha)) = ...
        real(n(activeLayer, find(wellenlange == alpha)))+1i*k; % substitute k of the active layer at one wavelength
                                                
    n = n(:, find(wellenlange == alpha));       % restrict optical constants to one wavelength
                                         
    wellenlange = alpha;
end
 
%--------------Standard transfer matrix simulation; for more information see 
%   Accurate Internal Quantum Efficiency Measurements in Organic and Other 
%   Thin Solar Cells, G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater., 22, 3293.
t(1)=0; 
t_cumsum=cumsum(t); 
x_pos=(1/2):1:sum(t); 
x_mat= sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1; 
R= wellenlange*0;
R_glass=abs((1-n(1,:))./(1+n(1,:))).^2;
E=zeros(length(x_pos),length(wellenlange));
for l = 1:length(wellenlange)
    S=I_mat(n(1,l),n(2,l)); 
    for matindex=2:(length(t)-1)
        S=S*L_mat(n(matindex,l),t(matindex),wellenlange(l))*I_mat(n(matindex,l),n(matindex+1,l)); % L_mat: layer matrix ~phase matrix! describing the propagation
    end
    R(l)=abs(S(2,1)/S(1,1))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
    T(l)=abs(2/(1+n(1,l)))/sqrt(1-R_glass(l)*R(l)); %Transmission of field through glass substrate Griffiths 9.85 + multiple reflection geometric series
 
    for devices = 2:length(t)
        xi=2*pi*n(devices,l)/wellenlange(l); 
        dj=t(devices);
        x_indices=find(x_mat == devices); 
        x=x_pos(x_indices)-t_cumsum(devices-1); 
        S_prime=I_mat(n(1,l),n(2,l));
        for matindex=3:devices
            S_prime=S_prime*L_mat(n(matindex-1,l),t(matindex-1),wellenlange(l))*I_mat(n(matindex-1,l),n(matindex,l));
        end
        S_doubleprime=eye(2); % eye(n) returns an n-by-n matrix with 1's at the diagonal and zeros everywhere else
        for matindex=devices:(length(t)-1)
            S_doubleprime=S_doubleprime*I_mat(n(matindex,l),n(matindex+1,l))*L_mat(n(matindex+1,l),t(matindex+1), wellenlange(l));
        end
        E(x_indices,l)=T(l)*(S_doubleprime(1,1)*exp(-1i*xi*(dj-x))+S_doubleprime(2,1)*exp(1i*xi*(dj-x))) ./(S_prime(1,1)*S_doubleprime(1,1)*exp(-1i*xi*dj)+S_prime(1,2)*S_doubleprime(2,1)*exp(1i*xi*dj));
    end
end
% Absorption coefficient in cm^-1
a = zeros(length(t),length(wellenlange)); % matrix with 6 rows and as many columns as lambda (451) because we start from 350 nm
for matindex=2:length(t)
    a(matindex,:)=4*pi*imag(n(matindex,:))./(wellenlange*1e-7); % n is a matrix with the optical constants in the order of the layers and for all wavelength
end
Pos = find(x_mat == activeLayer);
% position dependent absorption rate in the active layer as a function of wavelength
AbsRate_activeLayer = repmat(a(activeLayer,:).*real(n(activeLayer,:)),length(Pos),1).*(abs(E(Pos,:)).^2); % a(matindex,:) executes for all the matindexes
% position independent absorption of the active layer as a function of
% wavelength (Absorption_activeLayer can take values between 0 and 1)
Absorption_activeLayer = sum(AbsRate_activeLayer,1)*1e-7;

function I = I_mat(n1,n2)
r=(n1-n2)/(n1+n2);
t=2*n1/(n1+n2);
I=[1 r; r 1]/t;
 
function L = L_mat(n,d,lambda)
xi=2*pi*n/lambda;
L=[exp(-1i*xi*d) 0; 0 exp(1i*xi*d)];

function Spectrum = LoadFile(filename,wellenlange) % load EQE data
Spectrum = zeros(length(wellenlange),length(filename));
for index = 1:length(filename)
    [data,~] = xlsread(filename{index});
    wavelength = data(:,1);
    Data = data(:,2);
    Spectrum(:,index) = interp1(wavelength,Data,wellenlange,'linear');
    % checks if the specified wavelength range for the interpolation is larger than the experimental EQE spectral range or if one EQE spectrum has less wavelength range than nny other; if so error is produced 
    if any(isnan(Spectrum)) == 1
        throw(Loadfile)
    end
end
 
function ntotal = LoadRefrIndex(name,wellenlange, nk_filename) % load optical constants
% IndRefr comprises the data and IndRefr_names comprises the names in the first row
[IndRefr,IndRefr_names] = xlsread(nk_filename); 
% find the columns that belongs to the indexed material
idx = cellfun(@(x) isequal(x,strcat(name,'_wavelengths')) | isequal(x,strcat(name,'_n')) | isequal(x,strcat(name,'_k')),IndRefr_names(1,:));
% m contains the optical constants of one material and is NaN for empty rows
m = IndRefr(:,idx);
% discriminate whether n and k for the same material are given for different wavelengths
if size(m,2) == 3 % n&k use the same wavelength column, hence total number of columns is 3 for the indexed material
    m(any(isnan(m),2),:) = []; % empty space is deleted
    file_wavelengths = m(:,1);
    n = m(:,2);
    k = m(:,3);
    n_interp=interp1(file_wavelengths, n, wellenlange, 'linear','extrap'); % real part of complex n is interpolated
    k_interp=interp1(file_wavelengths, k, wellenlange, 'linear', 'extrap'); % imaginary part of complex n is interpolated
else % n&k have separate wavelengths columns, hence total number of columns is 4 for the indexed material
    n = m(:,2);
    k = m(:,4);
    k(any(isnan(k),2),:) = []; % deletes empty space if k or n have different number of rows
    n(any(isnan(n),2),:) = [];
    n_interp=interp1(m(1:length(n),1), n, wellenlange, 'linear','extrap'); % real part of complex n is interpolated
    k_interp=interp1(m(1:length(k),3), k, wellenlange, 'linear','extrap'); % imaginary part of complex n is interpolated
end
%Returns complex index of refraction of ONE material at all wavelengths
ntotal = n_interp+1i*k_interp;


