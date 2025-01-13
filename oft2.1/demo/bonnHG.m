% Plot for our THG experiment

% For copyright and contacts see readme.txt.

tic
clear; clc; close all;
addpath('..')

%% Configuration settings
materials = {'silica', 'air'};
argument = struct('type', 'DIA', 'min', 0.3, 'max', 0.6);
lambda = 1064;

modeTask1 = struct('type', 'HYBRID', 'lambda', lambda, 'nu', 1, 'maxmode', 1);
modeTask2 = struct('type', {{'HYBRID', 'te', 'tm'}}, 'lambda', lambda/2, 'nu', [0,2], 'maxmode', inf);
modeTask3 = struct('type', 'HYBRID', 'lambda', lambda/3, 'nu', [1,3], 'maxmode', inf);
infomode = false;

%% File to save/load modes
filename = 'bonnHG.mat';

%% Load or calculate modes
try
    fprintf('%s: Loading modes from a file...\n', mfilename);
    load(filename, 'MODES', '-mat');
catch ME
    if strcmpi(ME.identifier, 'MATLAB:load:couldNotReadFile')
        fprintf('%s: Couldn''t load modes from a file, building them...\n', mfilename);
        MODES1 = buildModes(argument, struct('materials', {materials}), modeTask1, infomode);
        MODES2 = buildModes(argument, struct('materials', {materials}), modeTask2, infomode);
        MODES3 = buildModes(argument, struct('materials', {materials}), modeTask3, infomode);
    else
        rethrow(ME);
    end;
    MODES = [MODES1; MODES2; MODES3];
    save(filename, 'MODES', '-mat');
end;

%% Show modes
fprintf('%s: Number of modes: %g\n', mfilename, length(MODES));
localFigure = figure; hold on;

showModes(MODES); drawnow;

%% Intersection points
fprintf('%s: Search for intersections...\n', mfilename);
PMPOINTS = phaseMatchingPoints(MODES, false, infomode);
i = 1;
while i<numel(PMPOINTS)+1
    if isnan(PMPOINTS(i).harmonic)
        PMPOINTS(i) = [];
    else
        i = i+1;
    end
end

showPMPoints(PMPOINTS); drawnow

tElapsed = toc;
fprintf('Time: %g\n', tElapsed);

overlap = [overlapIntegral(PMPOINTS([2,4,5])); ...
    overlapIntegral(PMPOINTS([2,4,5]),pi/2)];

fprintf('Overlap for THG: %g, %g, %g for %s, %s, %s\n', max(abs(overlap), [], 1), pmpDescription(PMPOINTS(2), 'short'), pmpDescription(PMPOINTS(4), 'short'), pmpDescription(PMPOINTS(5), 'short'));
