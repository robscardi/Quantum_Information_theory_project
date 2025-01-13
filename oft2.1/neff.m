function [N A] = neff(poiType, poi, par, fibreSpec, task, infomode, getWholeMode)
% Returns n_eff for particular mode and argument value

% Copyright: (cc-by) Fibres team, Meschede Group, Uni Bonn, 2009-2010
% http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com

if nargin == 0
    clear; clc; close all;
    
    fprintf('%s called without arguments---demo mode\n', upper(mfilename));
    
    poiType = 'dia'; % RAD or WVL
    poi = 0.22; % points of interest (RAD or WVL)
    par = 0.96; % LAMBDA, if poiType is RAD, or A, if poiType is WVL
    % harmonic = 1;
    fibreSpec = struct('materials', {{'silica', 'air'}}, 'coreCladdingRatio', 5.6/125);
    % mat_out = 'air';
    task = struct('modetype', 'HYBRID', 'modeindex', [1 1]);
    infomode = false;
    getWholeMode = ~true;
end


%% Prepare data
assert(correctModeOrder(task.modetype, task.modeindex(1)));

argmin = 0.9 * min(poi);
argmax = 1.1 * max(poi);

poi = sort(poi);

%% Find all modes
argument = struct('type', poiType, 'min', argmin, 'max', argmax, 'harmonic', 1);
%[STARTNEFF, STARTPAR] = FindModes(argument, par, materials, modetype, nu, infomode);
if ~isfield(task, 'maxmode')
    task.maxmode = Inf;
end
[STARTNEFF, STARTPAR] = findModes(argument, par, fibreSpec, task, infomode);
numberofmodes = length(STARTNEFF);

ModeVsLambda = strcmpi(poiType, 'WVL');

if numberofmodes > 200
    fprintf('%s: Too many modes (%g), tracing is likely to fail\n', upper(mfilename), numberofmodes);
    N = NaN;
    A = NaN;
    fprintf('%s: Argument: %g, n_eff-value: %g\n', upper(mfilename), A(1), N(1));
    return;
elseif numberofmodes == 0 && ~ModeVsLambda
    fprintf('%s: No modes found, adding fundamental mode explicitly\n', upper(mfilename));
    A = poi;
    N = RefrIndex(fibreSpec.materials{1}, par);
    fprintf('%s: Argument: %g, n_eff-value: %g\n', upper(mfilename), A(1), N(1));
    return;
elseif numberofmodes == 0 && ModeVsLambda
    fprintf('%s: No modes found. If you need just one wavelength, switch to the RAD mode\n', upper(mfilename));
    A = NaN;
    N = NaN;
    fprintf('%s: Argument: %g, n_eff-value: %g\n', upper(mfilename), A(1), N(1));
    return;
end

if infomode
    fprintf('%s: Number of modes found: %g\n', upper(mfilename), numberofmodes);
end

%% Choose the requested mode
if task.modeindex(2) > numberofmodes
    ME = MException('NEFF:ModeIndexHigherThanModeNumber', 'Not enough modes found: found %g modes, requested mode (%g, %g)\n', numberofmodes, modeindex(1), modeindex(2));
    throw(ME);
end;

StartNeff = STARTNEFF(task.modeindex(2));
StartArg = STARTPAR(task.modeindex(2));

%% Trace the requested mode
density = max(3, numberofmodes / 1.5);

%[mode] = TraceMode(StartArg, StartNeff, argument, par, materials, ...
%    modeindex, density, modetype, infomode);
mode = traceMode(StartArg, StartNeff, argument, par, ...
                fibreSpec, task, density, infomode);

% StartArg, StartNeff, argument, par, materials, ...
%     modeindex, density, modetype, infomode

%% Find the requested points in it
tolx = 1e-10;
options = optimset('TolX', tolx);

if length(mode.ARG) == 1
    fprintf('%s: No additional points found during tracing, using the starting value as the result.\n', upper(mfilename));
    
    mode.ARG = [mode.ARG poi];
    sort(mode.ARG);
    mode.NEFF = ones(length(mode.ARG)) * mode.NEFF(1);
else
    for i=1:length(poi)
        if ~any(mode.ARG == poi(i)) % this poi value is not found yet
            ElementBefore = find(poi(i) > mode.ARG, 1, 'last');
            if isempty(ElementBefore)
                % The mode starts later than poi(i)
                ElementBefore = 0;
                mode.ARG = [poi(i), mode.ARG(ElementBefore+1:end)];
                mode.NEFF = [NaN, mode.NEFF(ElementBefore+1:end)];
            elseif ElementBefore == length(mode.ARG)
                % the mode ended too early
                mode.ARG = [mode.ARG poi(i:end)]; 
                mode.NEFF = [mode.NEFF ones(1, length(poi) - i + 1) * NaN]; 
            else
                assert(mode.ARG(ElementBefore + 1) > poi(i));
                neffBefore = mode.NEFF(ElementBefore);
                neffAfter = mode.NEFF(ElementBefore + 1);
                
                if isnan(neffBefore) || isnan(neffAfter)
                    mode.ARG = [mode.ARG(1:ElementBefore), poi(i), mode.ARG(ElementBefore+1:end)];
                    mode.NEFF = [mode.NEFF(1:ElementBefore), NaN, mode.NEFF(ElementBefore+1:end)];
                else
                    if ModeVsLambda
                        d = par;
                        g = @(neff) fibreMode(d, neff, poi, fibreSpec, task);
                    else
                        lambda = par;
                        g = @(neff) fibreMode(poi, neff, lambda, fibreSpec, task);
                    end;
                    
                    try
                        [root, fval, fzeroResult] = fzero(g, [neffBefore neffAfter], options); %#ok<ASGLU>
                    catch ME1
                        ME2 = MException('NEFF:fzeroException', ...
                            'FZERO raised exception for Point of Interest %g (poi # %g)\n', ...
                            poi(i), i);
                        ME2 = addCause(ME2, ME1);
                        throw(ME2);
                    end;
                    if fzeroResult ~= 1
                        ME = MException('NEFF:fzeroFailed', ...
                            'FZERO failed with result %g for Point of Interest %g (poi # %g)\n', ...
                            fzeroResult, poi(i), i);
                        throw(ME);
                    else
                        mode.ARG = [mode.ARG(1:ElementBefore), poi(i), mode.ARG(ElementBefore+1:end)];
                        mode.NEFF = [mode.NEFF(1:ElementBefore), root, mode.NEFF(ElementBefore+1:end)];
                    end;
                end;
            end;
        end;
    end;
end;

if getWholeMode
    A = mode.ARG;
    if ModeVsLambda
        N = mode.NEFF;
    else
        N = mode.NEFF;
    end;
else
    IX = zeros(1, length(poi));
    for i=1:length(poi)
        IX(i) = find(mode.ARG == poi(i), 1, 'first');
    end;
    A = mode.ARG(IX);
    N = mode.NEFF(IX);
end;

if infomode
    for i=1:length(A)
        fprintf('%s: poiType: %s, argument: %g, n_eff-value: %g\n', upper(mfilename), poiType, A(i), N(i));
    end;
end;