function [E, H] = modeField(lambda, d, n_eff, fibreSpec, task, R, PHI) 
% Returns the field (E, H) value(s) for the two-layer modes (HE/EH, TE, TM)
% 
% MODEFIELD calculates field according to table 12-3 from Snyder,
% Love "Optical Waveguide Theory", page 250.
%
% TODO: Check: Light is assumed to be linearly polarized along x (phi = 0).
%
% R, PHI are equally-sized arrays or matrices, defining the points, in
% which field must be calculated.
% NU is the first modal index. The second modal index is defined by the (a,
% neff) pair.
% D is the core diameter. % TODO: check for three-layer structure
% The function returns 3 matrices: 
%    E - electric field components (3D matrix, made up of 3 stacked matrices. 
%    Each of them has the same size as R and PHI and contains one component
%    of the electric field (ER, EPHI, EZ respectively)
%    H - magnetic field components (3D matrix, made of 3 stacked matrices. 
%    Each of them has the same size as R and PHI and contains one component
%    of the magnetic field (HR, HPHI, HZ respectively)

% (by) Dan, Karapetyan, Pritzkau, Wiedemann, 2009
% http://agmeschede.iap.uni-bonn.de | kotya.karapetyan@gmail.com

lambda = lambda / 1000; % wavelength and geometry (d, R) should be in the same units

assert((size(R,1) == size(PHI,1)) && (size(R,2) == size(PHI,2)), ...
    '%s ERROR: size(R) ~= size(PHI)', mfilename);

assert(correctModeOrder(task.modetype, task.modeindex(1)), ...
    '%s: Invalide mode type (%s) or order (%g)\n', mfilename, ...
    task.modetype, task.modeindex(1));

a = d / 2; % fibre radius

nu = task.modeindex(1);
m = task.modeindex(2);

n_inn = refrIndex(fibreSpec.materials{1}, lambda);
n_out = refrIndex(fibreSpec.materials{2}, lambda);

%% Electromagnetic constants
eps0 = 8.8541878176 * 10^(-12); % permittivity of vacuum (SI)
mu0 = 1.2566 * 10^(-6); % permeability of vacuum (SI)

%% Propagation constant
k0 = 2*pi/lambda;
beta = n_eff * k0;

%% Auxiliary parameters
% Auxiliary parameters from Kien-2004.
h = sqrt(n_inn^2 * k0^2 - beta^2);
q = sqrt(beta^2 - n_out^2 * k0^2);
u = a * h;
w = a * q;
v = k0 * a * (n_inn^2 - n_out^2)^(1/2);

% Auxiliary parameters from Snyder/Love p. 250
b1 =   1 /(2*u) * (besselj(nu-1,u) / besselj(nu,u) - ...
    besselj(nu+1,u) / besselj(nu,u));
b2 = - 1 /(2*w) * (besselk(nu-1,w) / besselk(nu,w) + ...
    besselk(nu+1,w) / besselk(nu,w));
delta = (n_inn^2 - n_out^2) / 2 / n_inn^2; % p. 227, eq. 11-48
F1 = (u*w/v)^2 * (b1 + (1 - 2*delta) * b2) / nu;
F2 = (v/(u*w))^2 * nu / (b1 + b2);
a1 = (F2 - 1) / 2;
a2 = (F2 + 1) / 2;
a3 = (F1 - 1) / 2;
a4 = (F1 + 1) / 2;
a5 = (F1 - 1 + 2*delta) / 2;
a6 = (F1 + 1 - 2*delta) / 2;

if mod(m, 2) 
    % m is odd, e.g. for fundamental mode HE11 => HE modes 
    f = @(nu,phi) cos(nu*phi);
    g = @(nu,phi) -sin(nu*phi);
else
    % m is even => EH modes
    f = @(nu,phi) sin(nu*phi);
    g = @(nu,phi) cos(nu*phi);
end;

% Normalization of energy
%A = 1; % TODO: correct this, this is simplification
%g_in = A^2/(2*u);
%g_out = A^2*(besselj(v,h*a))^2 / (2*w*(besselk(v,q*a)^2));

%% Define functions for E and H 
% (all three vector components and both for inner and outer space)
if (strcmpi(task.modetype, 'HYBRID'))
    % Electric field inner (HE/EH)
    Ei_r   = @(r,phi) - (a1 * besselj(nu-1,u*r/a) + a2 * besselj(nu+1,u*r/a)) ./ besselj(nu,u) .* f(nu,phi);
    Ei_phi = @(r,phi) - (a1 * besselj(nu-1,u*r/a) - a2 * besselj(nu+1,u*r/a)) ./ besselj(nu,u) .* g(nu,phi);
    Ei_z   = @(r,phi) - (1i * u) / (a * beta) * besselj(nu,u*r/a) / besselj(nu,u) .* f(nu,phi);
    
    % Magnetic field inner (HE/EH)
    Hi_r   = @(r,phi)   sqrt(eps0/mu0) * (k0 * n_inn^2/beta) * (a3 * besselj(nu-1,u*r/a) - a4 * besselj(nu+1,u*r/a)) ./ besselj(nu,u) .* g(nu,phi);
    Hi_phi = @(r,phi) - sqrt(eps0/mu0) * (k0 * n_inn^2/beta) * (a3 * besselj(nu-1,u*r/a) + a4 * besselj(nu+1,u*r/a)) ./ besselj(nu,u) .* f(nu,phi);
    Hi_z   = @(r,phi) - sqrt(eps0/mu0) * (1i * u * F2) / (k0 * a) * besselj(nu,u*r/a) / besselj(nu,u) .* g(nu,phi);
    
    % Electric field outer (HE/EH)
    Eo_r   = @(r,phi) - u / w * (a1 * besselk(nu-1,w*r/a) - a2 * besselk(nu+1,w*r/a)) ./ besselk(nu,w) .* f(nu,phi);
    Eo_phi = @(r,phi) - u / w * (a1 * besselk(nu-1,w*r/a) + a2 * besselk(nu+1,w*r/a)) ./ besselk(nu,w) .* g(nu,phi);
    Eo_z   = @(r,phi) - (1i * u) / (a * beta) * besselk(nu,w*r/a) / besselk(nu,w) .* f(nu,phi);
    
    % Magnetic field outer (HE/EH)
    Ho_r   = @(r,phi)   sqrt(eps0/mu0) * (k0 * n_inn^2/beta) * (u / w) * (a5 * besselk(nu-1,w*r/a) + a6 * besselk(nu+1,w*r/a)) / besselk(nu,w) .* g(nu,phi);
    Ho_phi = @(r,phi) - sqrt(eps0/mu0) * (k0 * n_inn^2/beta) * (u / w) * (a5 * besselk(nu-1,w*r/a) - a6 * besselk(nu+1,w*r/a)) / besselk(nu,w) .* f(nu,phi);
    Ho_z   = @(r,phi) - sqrt(eps0/mu0) * (1i * u * F2) / (k0 * a) * besselk(nu,w*r/a) / besselk(nu,w) .* g(nu,phi);

elseif (strcmpi(modetype, 'TE'))
    % Electric field inner (TE)
    Ei_r   = @(r,phi)   0 .*r.*phi;
    Ei_phi = @(r,phi) - besselj(1,u*r/a) / besselj(1,u);
    Ei_z   = @(r,phi)   0 .*r.*phi;
    
    % Magnetic field inner (TE)
    Hi_r   = @(r,phi) sqrt(eps0/mu0) * (beta/k0) * besselj(1,u*r/a) / besselj(1,u);
    Hi_phi = @(r,phi) 0 .*r.*phi;
    Hi_z   = @(r,phi) sqrt(eps0/mu0) * (1i * u) / (k0 * a) * besselj(0,u*r/a) / besselj(1,u);
    
    % Electric field outer (TE)
    Eo_r   = @(r,phi)   0 .*r.*phi;
    Eo_phi = @(r,phi) - besselk(1,w*r/a) / besselk(1,w);
    Eo_z   = @(r,phi)   0 .*r.*phi;
    
    % Magnetic field outer (TE)
    Ho_r   = @(r,phi)   sqrt(eps0/mu0) * (beta/k0) * besselk(1,w*r/a) / besselk(1,w);
    Ho_phi = @(r,phi)   0 .*r.*phi;
    Ho_z   = @(r,phi) - sqrt(eps0/mu0) * (1i * w) / (k0 * a) * besselk(0,w*r/a) / besselk(1,w);

elseif (strcmpi(modetype, 'TM'))
    % Electric field inner (TM)
    Ei_r   = @(r,phi) besselj(1,u*r/a) / besselj(1,u);
    Ei_phi = @(r,phi) 0 .*r.*phi;
    Ei_z   = @(r,phi) 1i * u / (a * beta) * besselj(0,u*r/a) / besselj(1,u);
    
    % Magnetic field inner (TM)
    Hi_r   = @(r,phi) 0 .*r.*phi;
    Hi_phi = @(r,phi) sqrt(eps0/mu0) * (k0 * n_inn^2/beta) * besselj(1,u*r/a) / besselj(1,u);
    Hi_z   = @(r,phi) 0 .*r.*phi;
    
    % Electric field outer (TM)
    Eo_r   = @(r,phi) n_inn^2 / n_out^2 * besselk(1,w*r/a) / besselk(1,w);
    Eo_phi = @(r,phi) 0 .* r .*phi;
    Eo_z   = @(r,phi) -1i * n_inn^2 / n_out^2 * w / (a * beta) * besselk(0,w*r/a) / besselk(1,w); 
    
    % Magnetic field outer (TM)
    Ho_r   = @(r,phi) 0 .* r .*phi;
    Ho_phi = @(r,phi) sqrt(eps0/mu0) * (k0 * n_inn^2/beta) * besselk(1,w*r/a) / besselk(1,w);
    Ho_z   = @(r,phi) 0 .* r .* phi;

else
    throw(MException('FIELD:InvalidModeType'));
end;

%% Calculate fields
try
    if sum(sum(R > a)) == 0 % only R <= a present --> inner mode
        ER = Ei_r(R, PHI); 
        EPHI = Ei_phi(R,PHI);
        EZ = Ei_z(R,PHI);
        HR = Hi_r(R, PHI); 
        HPHI = Hi_phi(R,PHI);
        HZ = Hi_z(R,PHI);
    elseif sum(sum(R < a)) == 0 % only R >= a present --> outer mode
        ER = Eo_r(R, PHI);
        EPHI = Eo_phi(R, PHI);
        EZ = Eo_z(R, PHI);
        HR = Ho_r(R, PHI);
        HPHI = Ho_phi(R, PHI);
        HZ = Ho_z(R, PHI);
    else % both inner and outer space points present
        % Calculate for inner
        [row col] = find(R <= a);
        ER(row, col) = Ei_r(R(row, col), PHI(row, col)); 
        EPHI(row, col) = Ei_phi(R(row, col), PHI(row, col));
        EZ(row, col) = Ei_z(R(row, col), PHI(row, col));
        HR(row, col) = Hi_r(R(row, col), PHI(row, col)); 
        HPHI(row, col) = Hi_phi(R(row, col), PHI(row, col));
        HZ(row, col) = Hi_z(R(row, col), PHI(row, col));

        % Calculate for outer        
        [row col] = find(R > a);
        ER(row, col) = Eo_r(R(row, col), PHI(row, col));
        EPHI(row, col) = Eo_phi(R(row, col), PHI(row, col));
        EZ(row, col) = Eo_z(R(row, col), PHI(row, col));
        HR(row, col) = Ho_r(R(row, col), PHI(row, col));
        HPHI(row, col) = Ho_phi(R(row, col), PHI(row, col));
        HZ(row, col) = Ho_z(R(row, col), PHI(row, col));
    end;
    
    % Combine radial, azimuthal and longitudinal components of E and H
    % fields
    E = cat (3, ER, EPHI, EZ);
    H = cat (3, HR, HPHI, HZ);
    
catch ME 
    rethrow(ME);
end;


