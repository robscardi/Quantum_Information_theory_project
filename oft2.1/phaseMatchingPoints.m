function [PMPOINTS] = ...
    phaseMatchingPoints(MODES, sameWavelengthPoints, infomode)
% Finds all phase-matching points for the given modes
%
% If SAMEWAVELENGTHPOINTS is TRUE, will find intersection points also for
% the curves of the same wavelength.
%
% Uses the HARMONIC parameter of the modes to bring all modes up to their
% fundamental wavelength.

% Copyrtight (cc-by) Karapetyan, Pritzkau. Microfibres Team, AG Meschede,
% Uni Bonn.
% http://agmeschede.iap.uni-bonn.de | kotya.karapetyan@gmail.com


%% SETTINGS
if nargin == 0 % debug mode
    clear; clc; close all;
    
    mat_inn = 'silica';
    mat_out = 'air';

    argtype = 'RAD'; % either 'RAD' or 'WVL'
    argmin = 0.0001; % [um]
    argmax = 0.5; % [um]
    NU = [1]; %#ok<NBRAK>
    PAR = [2 0.5]; % [um], this is lambda, since argtype is RAD
    MODETYPE = {'HYBRID'};
    harmonic = 1;
    
    fprintf('%s: Infomode: Building modes...\n', mfilename);
    MODES = BuildModes(argmin, argmax, argtype, PAR, harmonic, ...
        mat_inn, mat_out, MODETYPE, NU, false);
   
    infomode = true;
    sameWavelengthPoints = false;
    
    PhaseMatchingPoints(MODES, sameWavelengthPoints, infomode);
    return;
end; 

if infomode
    fprintf('%s: Number of modes: %g\n', mfilename, length(MODES));
    oldFigure = get(0,'CurrentFigure');
    localFigure = figure; hold on;
    ShowModes(MODES);
end;

assert(nargin == 3);

roughIntersectionsCount = 0;

% Check argtype
argtype = MODES(1).argtype;
% Check that ARGTYPE is valid
assert(strcmpi(argtype, 'DIA') || strcmpi(argtype, 'WVL'), ...
    '%s: ERROR: Wrong mode argument type: %s\n', mfilename, upper(argtype));

ModeVsLambda = strcmpi(argtype, 'WVL');
for i=2:length(MODES)
    assert(ModeVsLambda == strcmpi(MODES(i).argtype, 'WVL'), ...
        'All modes should have equal argument type');
end;

PMPOINTS = [];

for i = 1:length(MODES)-1
    for j = i+1:length(MODES)
        MODE1 = MODES(i);
        MODE2 = MODES(j);
        if infomode
            figure(localFigure);
            mode1handle = plot(MODE1.ARG, MODE1.NEFF, '-', 'LineWidth', 2, 'Color', 'r');
            mode2handle = plot(MODE2.ARG, MODE2.NEFF, '-', 'LineWidth', 2, 'Color', 'r');
            fprintf('%s: (%g,%g): ', mfilename, i, j);
        end
        if ~sameWavelengthPoints % modes with same wavelengths should be excluded
            skip = false;
            if ModeVsLambda
                if MODE1.harmonic == MODE2.harmonic
                    skip = true;
                end;
            else
                if MODE1.par == MODE2.par
                    skip = true;
                end;
            end;
            if skip
                if infomode
                    fprintf('skipping modes with the same wavelengths\n');
                    try delete(mode1handle); delete(mode2handle); catch; end; %#ok<CTCH>
                end;
                continue;
            end;
        end;
        try
           [pmpoint] = phaseMatch(MODE1, MODE2, infomode);
           exitcode = pmpoint.quality;
        catch ME
            fprintf('%s: PHASEMATCH CALL FAILED: %s\n', upper(mfilename), ME.message);
            rethrow(ME);
        end;
        if abs(exitcode) == 1 % phase matching point found, fine or rough
            if exitcode == 1
                colour = 'k';
            else
                colour = 'r';
                roughIntersectionsCount = roughIntersectionsCount + 1;
                fprintf('%s: ROUGH RESULT\n', upper(mfilename));
            end;
            if infomode
                figure(localFigure);
                plot(pmpoint.arg, pmpoint.neff1, 'Marker', 'o', 'Color', colour, 'MarkerSize', 10);
                drawnow;
            end;
            PMPOINTS = [PMPOINTS; pmpoint]; %#ok<AGROW>
        end;
        if infomode; delete(mode1handle); delete(mode2handle); end;
    end;
end;

for i=2:length(PMPOINTS)
    for j=1:length(PMPOINTS)-1
       if PMPOINTS(i).arg < PMPOINTS(j).arg
           pmpoint = PMPOINTS(i);
           PMPOINTS(i) = PMPOINTS(j);
           PMPOINTS(j) = pmpoint;
       end;
    end;
end

if infomode
    fprintf('%s: Intersections: %g, of them rough: %g\n', ...
        mfilename, length(PMPOINTS), roughIntersectionsCount);
    
    try close(localFigure); catch; end %#ok<CTCH>
end;

try figure(oldFigure); catch; end; %#ok<CTCH>
