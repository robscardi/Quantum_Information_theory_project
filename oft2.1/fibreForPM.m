function PMPOINTS = fibreForPM(par, RANGE, materials, harmonic, displayModes, infomode)
% Returns possible phase-matching points for 2nd or 3rd harmonic generation
%
% EXAMPLES
%   fibreForPM(0.380, [840,1020], {'silica', 'air'}, 2, true, false) 
%   will search for 1st-to-2nd harmonics PM-points in a silica-in-air fibre
%   of 0.380 um diameter.
%   
%   fibreForPM(900, [0.3,0.5], {'silica', 'air'}, 3, true, false)
%   will search for 1st-to-3rd harmonics PM-points in a silica-in-air fibre
%   at 900 nm pump wavelength (THG at 300 nm)
%
%   fibreForPM(900) QUICK mode calculation for a given wavelength, SHG only
%
%   fibreForPM(0.5) QUICK mode calculation for a given diameter, SHG only

% Copyright: (cc-by) Subwavelength-diameter fibres team @ uni bonn, 2009-2010
% http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com

micronsThreshold = 200; % wavelength cannot be less than that in nanometers,
% diameter cannot be greater than that in microns.

if nargin == 1
    fprintf('QUICK MODE: ');
    if par < micronsThreshold % diameter given
        range = [830 1030];
        fprintf('will search for SHG phase-matching points between %g and %g nm, for silica-in-air.\n', range(1), range(2));
        fibreForPM(par, range, {'silica', 'air'}, 2, false, false);
    else % wavelength given
        range = [0.3 0.55];
        fprintf('will search for SHG phase-matching points between %g and %g um, for silica-in-air.\n', range(1), range(2));
        fibreForPM(par, range, {'silica', 'air'}, 2, false, false);
    end
    return
end

assert(nargin == 6, 'Invalid number of arguments');

MODETYPE = {'HYBRID'; 'TE'; 'TM'};
NU = 0:10;
% close all;

if displayModes
    localFigure = figure;
end
  
assert(numel(RANGE) == 2);
RANGE = sort(RANGE);

if par < micronsThreshold 
    %% PAR is diameter
    assert(sum(RANGE > micronsThreshold) == 2);
    fprintf('Diameter given (%g um), search for possible wavelengths. Building modes...', par);
    if RANGE(1) < 700
        RANGE = RANGE * harmonic;
        fprintf('\nWavelength range corrected to fundamental\n');
    end
    
    modeTask = struct('diameter', par, 'type', {MODETYPE}, 'nu', NU, 'maxmode', inf);

    % pump modes
    argument = struct('min', RANGE(1), 'max', RANGE(2), 'type', 'WVL', 'harmonic', 1);
    modesPump = buildModes(argument, struct('materials', {materials}), modeTask, false);
    
    % harmonic modes
    argument = struct('min', RANGE(1), 'max', RANGE(2), 'type', 'WVL', 'harmonic', harmonic);
    modesHarm = buildModes(argument, struct('materials', {materials}), modeTask, false);

elseif par > micronsThreshold % it is nanometers, so par is lambda
    %% PAR is lambda
    if par < 700
        fprintf('%s: WARNING: Wavelength of %g is probably the %g harmonic. Wavelength corrected to fundamental.\n', mfilename, par, harmonic);
        par = par * harmonic;
    end;
    assert(sum(RANGE < micronsThreshold) == 2,'ERROR: Diameter should be specified in microns');

    fprintf('Wavelength given (%g nm), search for possible diameters. Building modes...', par);
    modeTask = struct('lambda', par, 'type', {MODETYPE}, 'nu', NU, 'maxmode', inf);
    argument = struct('min', RANGE(1), 'max', RANGE(2), 'type', 'DIA');
	% materials = {'silica', 'air'};
    fibre = struct('materials', {materials});
    modesPump = buildModes(argument, fibre, modeTask, false);
    modeTask = struct('lambda', par / harmonic, 'type', {MODETYPE}, 'nu', NU, 'maxmode', inf);
    modesHarm = buildModes(argument, fibre, modeTask, false);
else
    error('Cannot guess working mode');
end;

if displayModes 
    figure(localFigure);
    showModes([modesPump; modesHarm]);
end

fprintf(' searching for phase-matching points...');
PMPOINTS = phaseMatchingPoints([modesHarm; modesPump], false, infomode);
fprintf(' done\n');
if displayModes 
    figure(localFigure);
    showPMPoints(PMPOINTS);
end
fprintf('\nThe following %g PM-points found\n', length(PMPOINTS));
for i = 1:length(PMPOINTS)
    fprintf('\nPM%g:\n', i);
    fprintf('%s\n', pmpDescription(PMPOINTS(i)));
end;
