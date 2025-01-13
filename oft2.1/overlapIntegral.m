function overlap = overlapIntegral(pmp, rotationAngle)
% Overlap integral for two modes at a phase-matching point 
%
% Note: Strictly speaking, phase-matching is not necessary for this
% function to work, it just needs PMP as a structure containing two modes
% and the point-of-interest. 
%
% See grubskysavchenko.m for examples.
%
% Copyright: (cc-by) Subwavelength-diameter fibres team @ Uni Bonn, 2009-2010 
% http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com
% This code was implemented in cooperation with Timothee Lee, ORC, Soton.

% TODO: Correct use of HARMONIC

% close all
if length(pmp) > 1
    for i = 1:length(pmp)
        if nargin == 1
            rotationAngle = 0;
        end
        overlap(i) = overlapIntegral(pmp(i), rotationAngle); %#ok<AGROW>
    end;
    return
end;

if nargin == 1 
    dPhi = 0;
else
    dPhi = rotationAngle;
end

%% Extract data from pmp structure
if strcmpi(pmp.argtype, 'wvl') % modes neff vs. lambda
    lambda1 = pmp.arg / pmp.mode1.harmonic;
%     harmonic1 = pmp.mode1.harmonic;
%     harmonic2 = pmp.mode2.harmonic;
    lambda2 = pmp.arg / pmp.mode2.harmonic;
    d = pmp.mode1.par;
    assert(d == pmp.mode2.par, '%s: diameters cannot be different\n', upper(mfilename));
else % modes neff vs. fibre diameter
    lambda1 = pmp.mode1.par;
    lambda2 = pmp.mode2.par;
%     if lambda1 > lambda2 % todo: Enable harmonics for this mode too. Make par real wavelength or PAR * HARMONIC real wavelength
%         harmonic1 = 1;
%         harmonic2 = round(lambda1 / lambda2 * 1e6) / 1e6; % integer harmonic if possible
%     else
%         harmonic1 = round(lambda2 / lambda1 * 1e6) / 1e6; % integer harmonic if possible
%         harmonic2 = 1;
%     end;
    
    d = pmp.arg;
end;

mat_inn = pmp.mode1.fibreSpec.materials{1};
mat_out = pmp.mode1.fibreSpec.materials{2};
assert(strcmpi(mat_inn, pmp.mode2.fibreSpec.materials{1}), '%s: Materials for intersecting modes cannot be different\n', upper(mfilename));
assert(strcmpi(mat_out, pmp.mode2.fibreSpec.materials{2}), '%s: Materials for intersecting modes cannot be different\n', upper(mfilename));

%% Check that this is THG, not SHG
if ~(pmp.mode2.harmonic == 3 && pmp.mode1.harmonic == 1)
    fprintf('%s: Only THG is supported\n', upper(mfilename));
    overlap = 0;
    return
end;

%% Define calculation grid
aperture = 9 * d;

Nr = [500, 1000]; % inn, out
Nphi = [32, 32]; % inn, out
dr_inn = d/2 / Nr(1);
dr_out = (aperture/2 - d/2) / Nr(2);

dphi_inn = 2 * pi / Nphi(1);
dphi_out = 2 * pi / Nphi(2);

F1 = struct('dr', [dr_inn dr_out], 'dphi', [dphi_inn dphi_out], ...
    'diam', [d aperture], 'E1', [], 'H1', [], 'E2', [], 'H2', []);
F2 = F1;
    
FG = fieldGrid(F1);

%% Calculate fields
% F1 is the first mode, F2 is the second one. E1 is the e-field inside the
% fibre, E2 is outside.
task1 = struct('modetype', pmp.mode1.modetype, 'modeindex', pmp.mode1.modeindex);
[F1.E1, F1.H1] = modeField(lambda1, d, pmp.neff1, pmp.mode1.fibreSpec, task1, FG.R1, FG.PHI1);
[F1.E2, F1.H2] = modeField(lambda1, d, pmp.neff1, pmp.mode1.fibreSpec, task1, FG.R2, FG.PHI2);
task2 = struct('modetype', pmp.mode2.modetype, 'modeindex', pmp.mode2.modeindex);
[F2.E1, F2.H1] = modeField(lambda2, d, pmp.neff2, pmp.mode2.fibreSpec, task2, FG.R1, FG.PHI1 + dPhi);
[F2.E2, F2.H2] = modeField(lambda2, d, pmp.neff2, pmp.mode2.fibreSpec, task2, FG.R2, FG.PHI2 + dPhi);

%% Make F1 the fundamental wave (bigger lambda)
if lambda1 < lambda2 % field 1 is acutally a harmonic of field 2
    % Shift fields
    F = F2;
    F2 = F1;
    F1 = F;
    clear F
end;

%% Normalize power
F1 = normalizePower(F1);
F2 = normalizePower(F2);

%% Z0 factor from Grubsky2005, p. 6799
eps0 = 8.8541878176e-12; % permittivity of vacuum (SI)
mu0  = 1.2566e-6; % permeability of vacuum (SI)
z0   = sqrt(mu0 / eps0);
% z0 = 1;

F1.E1 = F1.E1 / sqrt(z0);
F1.E2 = F1.E2 / sqrt(z0);
F2.E1 = F2.E1 / sqrt(z0);
F2.E2 = F2.E2 / sqrt(z0);

F1.H1 = F1.H1 * sqrt(z0);
F1.H2 = F1.H2 * sqrt(z0);
F2.H1 = F2.H1 * sqrt(z0);
F2.H2 = F2.H2 * sqrt(z0);

%% Plot the mode fields
% X = FG.R1 .* cos(FG.PHI1);
% Y = FG.R1 .* sin(FG.PHI1);
% DisplayFields(F1.E1(:,:,1), F1.E1(:,:,2), F1.E1(:,:,3), X, Y, d/2)
% DisplayFields(F2.E1(:,:,1), F2.E1(:,:,2), F2.E1(:,:,3), X, Y, d/2)


%% Overlap integral
J3_ScalarProduct_inn = sum(conj(F1.E1) .* F2.E1, 3) .* sum (conj(F1.E1) .* conj (F1.E1), 3) ;
J3_DiffScalarProduct_inn = J3_ScalarProduct_inn .* FG.ds1;
J3 = sum(sum(J3_DiffScalarProduct_inn));
% J3_ScalarProduct_out = dot (conj(F1.E2), F2.E2, 3) .* dot (conj(F1.E2), conj (F1.E2), 3) ;
% J3_DiffScalarProduct_out = J3_ScalarProduct_out .* ds_out;
% J3 = J3 + sum(sum(J3_DiffScalarProduct_out));

% fprintf('Overlap for modes %s and %s: %12.9g\n', ModeDescription(pmp.mode1), ModeDescription(pmp.mode2), J3);
overlap = J3;

%% Old Cristian's equations, calculated only inside the fibre (region 1)
% J1_ScalarProduct = 2 * (ModVector(F1.E1)).^4 + (ModVector(F1.E1.^2)).^2;
% J1_DiffScalarProduct = J1_ScalarProduct .* FG.ds1;
% J1 = sum(sum(J1_DiffScalarProduct)) / 3; %#ok<NASGU>
% 
% % Calculating J2 
% J2_ScalarProduct = ((ModVector(F1.E1)).^2) .* ((ModVector(F2.E1)).^2) + (abs(dot(F1.E1, F2.E1, 3))).^2 + (abs(dot(F1.E1, conj(F2.E1), 3))).^2;
% J2_DiffScalarProduct = J2_ScalarProduct .* FG.ds1;
% J2 = sum(sum(J2_DiffScalarProduct)) / 3; %#ok<NASGU>
% 
% % Calculating J5 
% J5_ScalarProduct = 2 * (ModVector(F2.E1)).^4 + (ModVector(F2.E1.^2)).^2;
% J5_DiffScalarProduct = J5_ScalarProduct .* FG.ds1;
% J5 = sum(sum(J5_DiffScalarProduct)) / 3; %#ok<NASGU>
% 
% function [A] = ModVector(E)
% 
% A = sqrt((abs(E(:,:,1))).^2 + (abs(E(:,:,2))).^2 + (abs(E(:,:,3))).^2);