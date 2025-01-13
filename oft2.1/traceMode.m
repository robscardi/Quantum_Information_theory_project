function [res] = traceMode(StartArg, StartNeff, argument, par, fibreSpec, ...
    task, density, infomode)
% Calculates the mode curve (n_eff vs d or lambda) starting from one point
%
% Finds the NEFF vs. ARG guided mode of a fibre starting from the
% point [STARTARG, STARTNEFF] which should belong to this mode and which
% determines which of the multiple possible roots of the modal equation
% (Tong et al. 2004. Eqs (3-5)) should be taken. MODEINDEX = [nu m] (mode
% indices). MAT_INN and MAT_OUT are the names of core and cladding (or
% material and surrounding for pulled fibre waist) as required by 
% FIBREINDEX 
%
% MODETYPE specifies the mode to be calculated: 'HYBRID', 'TE' or 'TM'
%
% See also BUILDMODES

% Copyright: (by) Karapetyan, Dan, Pritzkau, Wiedemann. 2008-2010. http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com

%% Tuning settings
ReduceStepFactor = 3; % how much the step should be decreased if needed, per cycle
n_shift_factor = 5; % the more  -- the less chance that step decrease will be needed, but the smaller approach to function zero
dampingweight = 10; % how much current step should damp returning to desired_step, from 0 to inf

pointcountmin = 50; % minimum amount of points

%% Input arguments check
assert(nargin == 8, 'Wrong number of input arguments');

argtype = argument.type;
assert(strcmpi(argtype, 'WVL') || strcmpi(argtype, 'DIA'), ...
    'Wrong mode parameter type %s', upper(argtype));

ModeVsLambda = strcmpi(argtype, 'WVL');

if strcmpi(argument.type, 'wvl')
    harmonic = argument.harmonic;
else
    harmonic = 1;
end

if ModeVsLambda
    StartArg = StartArg / harmonic;
    argument.min = argument.min / harmonic;
    argument.max = argument.max / harmonic;
end;

% assert(ismember(tracedirection, -1:1));

if infomode
    oldFigure = get(0,'CurrentFigure');
    localFigure = figure;
end;

%% Step settings
oddmode = mod(task.modeindex(2), 2) == 1;

desired_step = (argument.max - argument.min) / density;
if strcmpi(task.modetype, 'zhang') 
    % TODO: Make this smarter:
    % This is a workaround (2011-10-18), because without it the Zhang core
    % modes do not work properly. However, since they do not work anyway
    % (as of 2011-10-18), these lines might be removed. However, it would
    % be reasonable to change the desired_step in the core and cladding 
    % regions    
    desired_step = desired_step / 2;
end
% stepmin = desired_step / (ReduceStepFactor^20 + 1); % the smaller the slower and the finer approach to the end of the mode
stepmin = 1e-7; % step is always greater than dn, and dn < 1e-6 makes not much sense
if strcmpi(task.modetype, 'erdogan') 
    stepmin = 1e-4;
end

%% Refractive index vs. layer number
n = @(i, lambda) refrIndex(fibreSpec.materials{i}, lambda);

%% Set the inner and outer layer numbers
if sum(strcmpi(task.modetype, {'monerie', 'tsaohybrid', 'tsaote', 'tsaotm', 'erdogan', 'zhang'})) == 1
    innerLayer = 1;
    outerLayer = 3;
    assert(numel(fibreSpec.materials) == 3);
elseif sum(strcmpi(task.modetype, {'hybrid', 'te', 'tm'})) == 1
    switch numel(fibreSpec.materials)
        case 2
            innerLayer = 1;
            outerLayer = 2;
        case 3
            switch lower(task.region)
                case 'core'
                    innerLayer = 1;
                    outerLayer = 2;
                case 'cladding'
                    innerLayer = 2;
                    outerLayer = 3;
                otherwise
                    error('Wrong region spec');
            end
        otherwise
            error('Wrong number of materials');
    end
else
    error('Wrong mode type');
end 

%% n_max, n_min, F
if ModeVsLambda
    d = par;
    F = @(lambda, x) fibreMode(d, x, lambda, fibreSpec, task);
    n_max = @(x) n(innerLayer, x);
    n_min = @(x) n(outerLayer, x);
    % Check that n_max is above n_min more or less everywhere between
    % PARMIN and PARMAX
%     x = linspace(argmin, argmax, 100);
%     assert(sum(n_max(x) > n_min(x)) == length(x));
else % argtype == 'DIA'
    lambda = par;
    F = @(a, x) fibreMode(a, x, lambda, fibreSpec, task);
    % We don't really need a function here, but to make it compatible with
    % calls in WVL-mode:
    n_max = @(x) n(innerLayer, lambda);
    n_min = @(x) n(outerLayer, lambda);
    assert(n_max(0) > n_min(0));
end;

if infomode 
    plot(StartArg, StartNeff, 'ks', 'MarkerSize', 5);
end;

%% MODE TRACING
NEFF = [];
ARG = [];

nold = NaN;
argold = NaN;
step = 2 * stepmin;

% Start from the starting value and go left, then right
for direction = -1:2:1
    % Prepare for tracing
    %     if tracedirection ~= 0 && direction ~= tracedirection
    %         continue;
    %     end;
       
    nnew = StartNeff;
    argnew = StartArg;
    
    % Do tracing
    while step > stepmin % && ...
        %             abs(nold - n_min(argnew)) > 2 * eps && ...
        %             abs(argold - argmin) > 2 * eps && ...
        %             abs(argold - argmax) > 2 * eps
        if numel(NEFF) > 2
            if abs(NEFF(end-1) - NEFF(end)) < 1e-10
                break
                % this check is needed because step never reaches stepmin
                % as long as new points are successfully found
            end
        end
        
        if nnew ~= nold
            NEFF = [NEFF nnew]; %#ok<AGROW>
            ARG = [ARG argnew]; %#ok<AGROW>
            if infomode
                set(0,'CurrentFigure',localFigure);
                hold on;
                if ModeVsLambda
                    lbd = argnew;
                else
                    lbd = lambda;
                end;
                plot(argnew, nnew, '.', 'Color', colourVsLambda(lbd));
                try delete(circleHandle); catch; end; %#ok<CTCH>
                circleHandle = plot(argnew, nnew, 'ko', 'MarkerSize', 10);                
                drawnow;
            end;
            CalculateSlope = true;
            
            step = (dampingweight * step + desired_step) / (dampingweight + 1);
            
            nold = nnew;
            argold = argnew;
        end;
        
        if CalculateSlope
            if length(ARG) == 1
                % start with assumption of diagonal slope:
                slope = atan((n_max(argument.min) - n_min(argument.max)) / (argument.max - argument.min));
                if ModeVsLambda
                    slope = -slope;
                end;
            else
                % To continue, assume same slope as between last 2 points
                slope = atan((NEFF(end) - NEFF(end-1)) / (ARG(end) - ARG(end-1)));
            end;
            
            ssl = sin(slope); csl = cos(slope);
            CalculateSlope = false;
        end;
        
        darg = step * csl * direction; % arg-step
        argnew = argold + darg;
        % If arg is out of limits, set it to the limit:
        if argnew > argument.max
            step = step * (argument.max - argold) / darg;
            argnew = argument.max;
        end;
        if argnew < argument.min
            step = step * (argold - argument.min) / darg * direction;
            argnew = argument.min;
        end;
        
        dn = step * ssl * direction; % neff-step
        if abs(dn) > abs((n_max(argnew) - n_min(argnew)) / 500)
            step = step / ReduceStepFactor;
            continue
        end
            
        n_shift = abs(dn) / n_shift_factor; % for later use
        
        g = @(x) F(argnew, x); % modal function @ argnew
        
        % Odd modes tend to approach bottom smoothly, while even modes
        % "hit" it.
        if oddmode && abs(nold - n_min(argnew)) < (n_max(argnew) - n_min(argnew)) / 1000
            nLow = n_min(argnew) + 10*eps;
            nHigh = n_min(argnew) + (n_max(argnew) - n_min(argnew)) / 1000; % nold;
            try
                root = findRoot(g, [nLow, nHigh]);
                if sign(root - nold) == sign(slope * direction)
                    % Check that the function goes on to increase (or decrease)
                    nnew = root;
                    if abs(nnew - nold) < 1e-9
                        break;
                    end;
                    continue;
                else
                    step = step / ReduceStepFactor;
                end;
            catch ME %#ok<NASGU>
                step = step / ReduceStepFactor;
                continue;
            end;
        else
            nHigh = nold + dn + n_shift; % two guess-borders, one a little higer than the guess point,
            nLow = nold + dn - n_shift; % another a little lower
            nLow = max(n_min(argnew) + 10*eps, nLow); % in case nLow < n_min
            if nHigh <= nLow
                step = step / ReduceStepFactor;
                continue;
            end;
        end;
        
        try
            root = findRoot(g, [nLow, nHigh]);
            if sign(root - nold) == sign(slope * direction)
                % Check that the function goes on to increase (or decrease)
                nnew = root;
                continue;
            end;
        catch %#ok<CTCH>
        end;
        
        G = g([nLow nHigh]);
        if ~oddmode
            if nHigh < n_min(argnew)+1e-6
                if sign(g(n_min(argnew)+eps)) == sign(g(nHigh))
                    break;
                end;
            end;
        end;
        
        if ~isreal(G)
            step = step / ReduceStepFactor;
            continue;
        end;
        func_at_nLow = G(1); sL = sign(func_at_nLow);
        func_at_nHigh = G(2); sH = sign(func_at_nHigh);
        
        % If margins are not at two sides from function zero-point, move
        % them until it is the case.
        % This should only be done, if not horizontal approach, because in
        % that case there is no space to play, so we just stop.
        decreasestep = false;
        if (step > stepmin) && (sL == sH)
            % In which direction should we search for zero?
            if abs(func_at_nHigh) > abs(func_at_nLow)
                dir = -1; % go towards nLow
            elseif abs(func_at_nHigh) < abs(func_at_nLow)
                dir = 1; % go towards nHigh
            else
                % func_at_nHigh == func_at_nLow, we have no idea, in which direction to go
                step = step / ReduceStepFactor;
                continue;
            end
            
            % Search for change of sign
            decreasestep = false;
            while sL == sH % zero is not between nLow and nHigh
                assert(abs(dir) == 1); % check that dir is +1 or -1
                switch dir
                    case 1 % going towards nHigh
                        if abs(func_at_nHigh) > abs(func_at_nLow) %...but should go to nLow
                            % => local minimum, get out of here!
                            decreasestep = true;
                            break;
                        end;
                        % Move margins higher
                        nLow = nHigh;
                        func_at_nLow = func_at_nHigh;
                        sL = sH;
                        nHigh = nHigh + n_shift;
                        if nHigh > n_max(argnew)
                            nHigh = n_max(argnew);
                            break;
                        end;
                        func_at_nHigh = g(nHigh);
                        sH = sign(func_at_nHigh);
                    case -1 % going towards nLow
                        if abs(func_at_nHigh) < abs(func_at_nLow) %...and should go to nHigh
                            % => local minimum, get out!
                            decreasestep = true;
                            break;
                        end;
                        % Move margins lower
                        nHigh = nLow;
                        func_at_nHigh = func_at_nLow;
                        sH = sL;
                        nLow = nLow - n_shift;
                        if nLow < n_min(argnew)
                            nLow = n_min(argnew);
                            break;
                        end;
                        func_at_nLow = g(nLow);
                        sL = sign(func_at_nLow);
                end % case dir
                if infomode
                    assert(isreal(func_at_nLow) && isreal(func_at_nHigh));
                end;
                if nLow > nold + 1e-6;
                    break;
                end;
            end % while: Search for change of sign
        end
        
        if decreasestep
            step = step / ReduceStepFactor;
            continue;
        end;
        
        % Now the guess values are at two sides from zero. Try fzero
        % between them.
        try
            root = findRoot(g, [nLow nHigh]);
            if sign(root - nold) == sign(slope * direction)
                % Check that the function goes on to increase (or decrease)
                nnew = root;
            else
                step = step / ReduceStepFactor;
            end;
        catch ME  %#ok<NASGU>
            step = step / ReduceStepFactor;
            continue;
        end
    end; % while
    [ARG, XI] = sort(ARG);
    NEFF = NEFF(XI);
    
    nold = StartNeff;
    argold = StartArg;
    step = desired_step;
    
    CalculateSlope = true; 
    decreasestep = false; 
end % direction

[ARG, XI] = sort(ARG);
NEFF = NEFF(XI);

if infomode
    if length(ARG) < 2
        error('%s: LESS THAN 2 POINTS!\n', upper(mfilename));
    end;
    set(0,'CurrentFigure',localFigure);
    axis([min(ARG) max(ARG) min(NEFF) max(NEFF)]);
end;

res = struct(...
    'NEFF', NEFF, ...
    'ARG', ARG, ...
    'argtype', argtype, ...
    'harmonic', harmonic, ...
    'par', par, ...
    'fibreSpec', fibreSpec, ...
    'modeindex', task.modeindex, ...
    'modetype', task.modetype);

if isfield(task, 'region')
    res.foundin = lower(task.region);
else
    res.foundin = [];
end

if infomode
    try 
        delete(circleHandle); 
    catch ME; 
        if ~strcmpi(ME.identifier, 'MATLAB:UndefinedFunction'); 
            rethrow(ME)
        end
    end;
end

if infomode
    close(localFigure);
    try figure(oldFigure); catch; end;  %#ok<CTCH>
end;

%% Check if sufficient amount of points found
if length(res.NEFF) < 2
    fprintf('%s WARNING: too little amount of points in the mode\n', mfilename);
else
    % amount of points in mode is too little, let's increase
    while length(res.NEFF) < pointcountmin
        m = res;
        try
            res = addPointsToMode(m, false, infomode);
        catch ME
            sprintf('%s: WARNING! Adding points fails with message: %s\n', mfilename, ME.message)
            break
        end
    end;
end;

res.ARG = res.ARG * harmonic;

