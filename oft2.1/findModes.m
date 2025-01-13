function [NEFF, ARG, density] = findModes(argument, par, fibreSpec, task, infomode)
% Returns starting points n_eff(arg) for all modes found according to task
%
% Tries to find all guided modes of order NU within [ARGMIN; ARGMAX] for
% given materials. Variable argument type (wavelength or radius) is
% specified by ARGTYPE. C is either the wavelength (if argument is radius)
% or radius (if argument is wavelength). 
%
% A mode is a line NEFF vs. ARG, along which the FIBREMODE function
% equals zero. 
%
% DENSITY is the number of change-of-sign points in the region of interest.
% The modified function (2011-11) does not return all found modes, but
% only  the requested amount, to speed up. Therefore the number of modes
% available remains unknown outside the function. DENSITY is later used in
% buildModes.m.
%
% (by) Karapetyan, Pritzkau, Wiedemann, 2008-2009.
% http://agmeschede.iap.uni-bonn.de | kotya.karapetyan@gmail.com


% TODO: decide on use of HARMONIC. The idea is to always specify the
% correct data for wavelength and only use HARMONIC to show the modes (i.e.
% mode between 400 and 500 nm can be shown between 800 and 1000) and find
% intersections 

NEFF = [];
ARG = [];
density = [];

if ~correctModeOrder(task.modetype, task.modeindex(1))
    return
end

%% Input arguments check
assert(nargin == 5, 'Invalid number of input arguments');

assert(strcmpi(argument.type, 'WVL') || strcmpi(argument.type, 'DIA'), ...
    'Wrong mode argument type %s', upper(argument.type));

assert(argument.min ~= argument.max);

%% Tuning settings
tolx = 1e-10;
tolfun = 1e-3;

%% Work
ModeVsLambda = strcmpi(argument.type, 'WVL');

if strcmpi(argument.type, 'wvl')
    harmonic = argument.harmonic;
else
    harmonic = 1;
end

materials = fibreSpec.materials;

n = @(i, lambda) refrIndex(materials{i}, lambda);
% k0 = @(lambda) 2 * pi ./ lambda;

assert(task.modeindex(1) >= 0);

% Define outerLayer
if numel(fibreSpec.materials) == 2 || strcmpi(task.region, 'core')
    outerLayer = 2;
elseif strcmpi(task.region, 'cladding')
    outerLayer = 3;
else
    error('Cannot set outerLayer\n');
end

if ModeVsLambda
    d = par; 
    F = @(arg, neff) fibreMode(d, neff, arg, fibreSpec, task);
    n_min = n(outerLayer, argument.max);
    n_max = n(outerLayer - 1, argument.min);
else % argtype == 'DIA'
    lambda = par;
    F = @(arg, neff) fibreMode(arg, neff, lambda, fibreSpec, task);
    n_min = n(outerLayer, lambda);
    n_max = n(outerLayer - 1, lambda);
end;

if ModeVsLambda
    % Diagonal line goes from bottom left to top right
%     % beta vs. lambda
%     arg_vs_beta = @(beta) argmin + (beta - betamin) * (argmax - argmin) / (betamax - betamin);
    % n_eff vs. lambda
    arg_vs_neff = @(n) argument.min/harmonic + (n - n_min) * (argument.max/harmonic - argument.min/harmonic) / (n_max - n_min);
else
    % Diagonal line goes from top left to bottom right
    arg_vs_neff = @(n) argument.max - (n - n_min) * (argument.max - argument.min) / (n_max - n_min);
end;

% if ModeVsLambda
%     beta_vs_neff = @(n) k0(arg_vs_neff(n)) .* n;
% else
%     beta_vs_neff = @(n) 2*pi/lambda .* n;
% end;

% Plot the line from which the modes are traced
if infomode
    line(arg_vs_neff([n_min, n_max]), [n_min, n_max], 'Color', 'k');
    drawnow
end;

% Function along the diagonal line
g = @(neff) F(arg_vs_neff(neff), neff); % n_eff vs lambda

% Function tabulation
N = 1e4;
doCalculation = true;
while doCalculation
    %% First, find the change-of-sign points (possible root positions)
    step = -(n_max - n_min) / N; % n_eff vs lambda
    x = n_max: step: n_min; % n_eff vs. lambda
    y = g(x);
    
    % consider only real points 
    if ~isreal(y)
%         XIpos = find(angle(y) == 0);
%         XIneg = find(angle(y) == pi);
%         XI = [XIneg, XIpos];
%         XI = sort(XI);
%         y = y(XI);
%         x = x(XI);
        ix = sign(y) == 1 | sign(y) == -1;
        x = x(ix);
        y = y(ix);
    end; 
    
    % Function normalization
    m = max(abs(y));
    gnorm = @(x) g(x) / m;
    
    if infomode
        % Remember the current figure handle
        oldFigure = get(0,'CurrentFigure');
        % Plot function for ?graphical solution?
        localFigure = figure; hold on;
        plot(x, y, 'b.-');
        ylim([-0.1 .1]);
        if ModeVsLambda
            title(sprintf('Diam: %g', par));
        else
            title(sprintf('Wvl: %g', par));
        end;
        xlabel('n_{neff}'); % n_eff vs lambda
        ylabel('F (mode function)');
        drawnow
        grid
    end;
    
    % Search for places where function changes sign
    signchanged = abs(sign(y(2:end)) - sign(y(1:end-1)));
    
    % Create array of values of X where function changes sign
    L = x(signchanged > 0);
    clear('signchanged');
    if numel(L) > N / 200
        if infomode
            fprintf('%s: Not enough points (%g) for so many change-of-sign points (%g). Increasing.\n', ...
                mfilename, N, numel(L))
        end
        N = N * 10;
        if infomode
            close(localFigure);
            try figure(oldFigure); catch; end; %#ok<CTCH>
        end;
        continue % re-calculate L
    else
        doCalculation = false;
    end

    %% Now find actual roots
    % Options for FZERO
    options = optimset('TolX', tolx, 'TolFun', tolfun);
    
    % Search for zeros
    ROOTS = [];
    for i=1:length(L)
        assert(sign(gnorm(L(i))) ~= sign(gnorm(L(i)+step)))
        if infomode
            h1 = line([L(i) L(i)], [-0.05 0.05], 'Color','red','LineStyle','--');
            h2 = line([L(i)+step L(i)+step], [-0.05 0.05], 'Color','red','LineStyle','--');
            drawnow
        end;
        try
            [root, fval, result] = fzero(gnorm, [L(i) L(i) + step], options); %#ok<ASGLU>
        catch ME  %#ok<NASGU>
            if infomode
                delete(h1, h2)
            end
            continue;
        end;
        if result == 1 % if result ~= 1, it's not a real root
            ROOTS = [ROOTS root]; %#ok<AGROW>
        end;
        if infomode
            delete(h1, h2)
        end
        
        % If enough roots found, stop
        if numel(ROOTS) == task.maxmode
            break
        end
    end;
    
    if infomode
        fprintf('%s: Number of modes (%s) found: %g\n', mfilename, upper(task.modetype), length(ROOTS));
        
        % Show "graphical solution"
        for i = 1:length(ROOTS);
            line([ROOTS(i) ROOTS(i)], [-0.05 0.05], 'Color','red');
        end;
        
        grid on;
        drawnow;
        
        fprintf('%s: Time: %g s\n', mfilename, toc);
    end;
    
    NEFF = ROOTS; % n_eff vs lambda
end
if infomode
    fprintf('%s: Number of points used: %g\n', mfilename, N);
end
ARG = arg_vs_neff(ROOTS); % n_eff vs lambda

[NEFF, XI] = sort(NEFF, 'descend');
ARG = ARG(XI);

if ModeVsLambda && harmonic ~= 1
    ARG = ARG * harmonic;
end;

if infomode
    close(localFigure);
    try figure(oldFigure); catch; end; %#ok<CTCH>
end;

density = numel(L);
