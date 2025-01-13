clc
clear
close all

tStart = tic;

% Note: diameter is specified in micrometers, wavelength in nanometers

%% Specify the fibre parameters
% Fibre materials (core, cladding)
materials = {'sm800core'; 'silica'};

% Fibre structure description
fibre = struct(...
    'materials', {materials});

%% Create the task for dispersion curves calculation

% Argument for dispersion curves calculation
argument = struct(...
    'type', 'wvl',... % calculate vs. wavelength
    'harmonic', 1,... % required 
    'min', 600,... % calculate from
    'max', 2000); % calculate to

% Specify which modes to search for
task = struct(...
    'nu', [0 1],... % first modal index
    'type', {{'hybrid', 'te'}},... % mode types
    'maxmode', 1,... % how many modes of each type and NU to calculate
    'diameter', 8.2);%,... % parameter, structure diameter, if argument is wavelength
    %'region', 'cladding');

%% Find guided modes and calculate the dispersion curves
addpath('..')
fprintf('%s', 'Calculating... ');
modes = buildModes(argument, fibre, task, false);    
fprintf('%s\n', 'Done');

%% Print which modes have been calculated
fprintf('Modes found:\n');
for i = 1:numel(modes)
   fprintf('%s\n', modeDescription(modes(i), false));
end

%% Display calculated dispersion curves
h = showModes(modes, 'Modal dispersion in a fibre, modes vs. wavelength');
set(h(1), 'Color', 'b')
set(h(2), 'Color', 'r')
legend(modeDescription(modes(1), false, false), modeDescription(modes(2), false, false));

%% Convert modes to vs. V-parameter and show
modesV = modeVsV(modes);
    
figure
h = showModes(modesV, 'Modal dispersion in a fibre, modes vs. V-parameter');
set(h(1), 'Color', 'b')
set(h(2), 'Color', 'r')
legend(modeDescription(modes(1), false, false), modeDescription(modes(2), false, false));

%% Show refractive index lines
% Core
x = linspace(argument.min, argument.max, 100);
y = refrIndex(materials{1}, x);
h = line(vParameter(task.diameter, x, materials{1}, materials{2}), y);
set(h, 'Color', 'k')
set(h, 'LineStyle', '--')

% Cladding
x = linspace(argument.min, argument.max, 100);
y = refrIndex(materials{2}, x);
h = line(vParameter(task.diameter, x, materials{1}, materials{2}), y);
set(h, 'Color', 'k')
set(h, 'LineStyle', '--')

fprintf('Elapsed time: %g s\n', toc(tStart))
