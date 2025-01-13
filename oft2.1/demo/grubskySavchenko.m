%% Optical fibre toolbox tutorial: Grubsky-Savchenko
% Recreates the plot from Grubsky, Savchenko, 2005.

% For copyright and contacts see readme.txt.

clear; close all; % clc;
tStart = tic;
addpath('..')

%% Configuration settings
materials = {'silica', 'air'}; 
fibre.materials = materials;
argument = struct('type', 'DIA', 'min', 0.3, 'max', 3);
modeTask1 = struct('type', {'HYBRID'}, 'lambda', 1064, 'nu', 1, 'maxmode', 1);
modeTask3 = struct('type', {'HYBRID'}, 'lambda', 1064/3, 'nu', 1, 'maxmode', inf);
infomode = false;

%% File to save/load modes
filename = 'grubskysavchenkomodes.mat';

%% Load or calculate modes
try
    fprintf('%s: Loading modes from a file...\n', mfilename);
    load(filename, 'MODES', '-mat');
catch ME
    if strcmpi(ME.identifier, 'MATLAB:load:couldNotReadFile')
        fprintf('%s: Couldn''t load modes from a file, building them...\n', mfilename);
        MODES1 = buildModes(argument, fibre, modeTask1, infomode);
        MODES3 = buildModes(argument, fibre, modeTask3, infomode);
    else
        rethrow(ME);
    end;
    MODES = [MODES1; MODES3];
    save(filename, 'MODES', '-mat');
end;

%% Show modes
fprintf('%s: Number of modes: %g\n', mfilename, length(MODES));
localFigure = figure; hold on;

showModes(MODES); drawnow;

%% Intersection points
fprintf('%s: Search for intersections...\n', mfilename);
PMPOINTS = phaseMatchingPoints(MODES, false, infomode);

showPMPoints(PMPOINTS); drawnow

roughIntersectionsCount = 0;
for i=1:length(PMPOINTS)
    if PMPOINTS(i).quality == -1
        roughIntersectionsCount = roughIntersectionsCount + 1;
    end;
end;

fprintf('%s: SEARCH DONE. Intersections: %g, of them rough: %g\n', ...
     mfilename, length(PMPOINTS), roughIntersectionsCount);

fprintf('Overlap integral from Grubsky-Savchenko for all found intersection points:\n');
overlapIntegral(PMPOINTS)
fprintf('Rotated by 90 degrees:\n');
overlapIntegral(PMPOINTS, pi/2)
 
fprintf('Time: %g\n', toc(tStart));

% return

% %% Check that the phase-matching points are unchanged
% check = pmCheckValue(PMPOINTS);
% % if check == 4116359 % silica/air, 0.3...3, HYBRID, 1.064, NU=[1 3]
% if check == 2274280 % silica/air, 0.3...3, HYBRID, 1.064, NU=1
% % if check == 2297044 % silica/air, 0.3...3, HYBRID, 1.064, NU=1, w/o overlap
%     fprintf('\n%s: PM result verified\n', mfilename);
% else
%     fprintf('\n%s: PM RESULT NOT VERIFIED\n', mfilename);
% end;
% 
% %% Check that the modes are unchanged
% check = modeCheckValue(MODES);
% % if check == 4050492 % silica/air, 0.3...3, HYBRID, 1.064, NU=[1 3]
% if check == 4108890 % silica/air, 0.3...3, HYBRID, 1.064, NU=1
%     fprintf('\n%s: Modes result verified\n', mfilename);
% else
%     fprintf('\n%s: MODES RESULT NOT VERIFIED\n', mfilename);
% end;
% 
