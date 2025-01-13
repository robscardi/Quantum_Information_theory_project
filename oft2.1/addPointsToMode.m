function res = addPointsToMode(mode, range, infomode, localFigure)
% Adds points to the mode curve, trying to double their density 
% in the RANGE of argument. RANGE limits the argument; if missing or FALSE,
% the full mode argument span is used as range.
%
% The function inserts a new point into the middle between the points,
% between which the distance is larger than the average inter-point
% distance divided by 2. Thus, the new points are inserted in the
% low-density places.

if numel(mode) > 1
    error('%s does not support multiple mode processing yet.\n', mfilename)
end

NEFF = mode.NEFF;
ARG = mode.ARG;
ModeVsLambda = strcmpi(mode.argtype(1), 'w');

if ~exist('range', 'var') || (range(1) == false && range(end) == false)
    range = [ARG(1) ARG(end)];
end
first = find(ARG<=range(1),1,'last');
if isempty(first)
    first = 1;
end
last = find(ARG>=range(2),1,'first');
if isempty(last)
    last = numel(NEFF);
end

if ~exist('infomode', 'var')
    infomode = false;
elseif infomode
    if ~exist('localFigure', 'var')
        oldFigure = get(0,'CurrentFigure');
        localFigure = figure;
        plot(ARG, NEFF, '.')
        ylabel('NEFF')
        if ModeVsLambda
            xlabel('Wavelength, nm')
        else
            xlabel('Diameter, um')
        end
        hold on
        xlim(range)
        ylim([min([NEFF(first), NEFF(last)]), max([NEFF(first), NEFF(last)])])
    end
end

n = @(i, lambda) refrIndex(mode.fibreSpec.materials{i}, lambda);
oddmode = mod(mode.modeindex(2), 2) == 1;

% Set the inner and outer layer numbers
if strmatch(lower(mode.modetype), strvcat('monerie', 'tsaohybrid', 'tsaote', 'tsaotm', 'erdogan', 'zhang'))
    innerLayer = 1;
    outerLayer = 3;
    assert(numel(mode.fibreSpec.materials) == 3);
elseif strmatch(lower(mode.modetype), strvcat('hybrid', 'te', 'tm'))
    switch numel(mode.fibreSpec.materials)
        case 2
            innerLayer = 1;
            outerLayer = 2;
        case 3
            switch lower(mode.foundin)
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

localTask = struct('modetype', mode.modetype, 'region', mode.foundin, 'modeindex', mode.modeindex);
if ModeVsLambda
    d = mode.par;
s    n_max = @(x) n(innerLayer, x);
    n_min = @(x) n(outerLayer, x);
    deltaN = @(x) n_max(x) - n_min(x);
    % Check that n_max is above n_min more or less everywhere between
    % PARMIN and PARMAX
%     x = linspace(argmin, argmax, 100);
%     assert(sum(n_max(x) > n_min(x)) == length(x));
else % argtype == 'DIA'
    lambda = mode.par;
    F = @(a, x) fibreMode(a, x, lambda, mode.fibreSpec, localTask);
    % We don't really need a function here, but to make it compatible with
    % calls in WVL-mode:
    n_max = @(x) n(innerLayer, lambda);
    n_min = @(x) n(outerLayer, lambda);
    deltaN = @(x) n_max(x) - n_min(x);
    assert(n_max(0) > n_min(0));
end;

DN = NEFF(2:end) - NEFF(1:end-1);
DA = ARG(2:end) - ARG(1:end-1);
DIST2 = DN.^2 + DA.^2;
DIST = sqrt(DIST2);
threshold = mean(DIST(first:last-1)) / 2;

NEWARGS = ones(1, length(ARG)-1) * NaN;
NEWNEFF = ones(1, length(ARG)-1) * NaN;
for i=first:last-1
    if infomode
        figure(localFigure); % bring focus to local figure in case it was switched to smth. else
        try delete(circleHandle); catch; end; %#ok<CTCH>
        try delete(circleHandle1); delete(circleHandle2); catch; end; %#ok<CTCH>
        circleHandle1 = plot(ARG(i), NEFF(i), 'ko', 'MarkerSize', 10);
        circleHandle2 = plot(ARG(i+1), NEFF(i+1), 'ko', 'MarkerSize', 10);
        drawnow;
    end;
    
    argnew = (ARG(i)+ARG(i+1)) / 2;
    
    if DN(i) < deltaN(argnew)/20 && DIST(i) < threshold 
        continue;
    end;
    
    if infomode
        figure(localFigure); % bring focus to local figure in case it was switched to smth. else
        try  delete(lineVertHandle); catch; end; %#ok<CTCH>
        lineVertHandle = line([argnew argnew], [n_min(argnew) n_max(argnew)], 'Color','red');
        drawnow;
    end;
    
    g = @(x) F(argnew, x);
    try
        if ~oddmode || NEFF(i) > n_min(argnew) + (n_max(argnew) - n_min(argnew)) / 10000
            %b = 1;
            limits = NEFF(i:i+1);
        else
            %b = 2;
            limits = [n_min(argnew)+eps, max(NEFF(i:i+1))];
        end;
        if infomode
            figure(localFigure); % bring focus to local figure in case it was switched to smth. else
            lineTopHandle = line([min(ARG) max(ARG)], [limits(1) limits(1)], 'Color','red');
            lineBottomHandle = line([min(ARG) max(ARG)], [limits(2) limits(2)], 'Color','red');
            delete(circleHandle1); delete(circleHandle2);
            drawnow;
        end;
        newneff = findRoot(g, limits);
    catch ME
        if i == 1 || isnan(NEFF(i))
            NEFF(i) = NaN;
            ARG(i) = NaN;
        else
            if strncmpi(ME.message, 'FINDROOT', 8) % FIND ROOT failed
                error('%s: Failed to add point at %g\n', mfilename, argnew)
            else
                rethrow(ME)
            end
        end
%             %                     rethrow(ME);
%             stopAtArg = ARG(i);
%             if infomode
%                 fprintf('%s: FIND ROOT FAILED. REMOVING ALL THE POINTS AFTER THIS ONE', upper(mfilename));
%                 hold off
%                 plot(ARG(1:i), NEFF(1:i), '.', 'Color', [0.5 0 1]);
%                 hold on
%                 plot(NEWARGS(1:end), NEWNEFF(1:end), '.', 'Color', [0.5 0 1]);
%                 plot(ARG(i), NEFF(i), 'rs', 'MarkerSize', 15);
%             end;
%             break;
%         end;
        continue;
    end;
    NEWNEFF(i) = newneff;
    NEWARGS(i) = argnew;
    if infomode
        figure(localFigure); % bring focus to local figure in case it was switched to smth. else
        plot(argnew, newneff, '.', 'Color', [0.5 0 1]);
        circleHandle = plot(argnew, newneff, 'ko', 'MarkerSize', 10);
        delete(lineTopHandle); delete(lineBottomHandle);
        delete(lineVertHandle);
        drawnow;
        hold on;
    end;
    
end;
[ARG, IX] = sort([ARG NEWARGS]);
NEFF = [NEFF NEWNEFF];
NEFF = NEFF(IX);
ARG(isnan(ARG)) = [];
NEFF(isnan(NEFF)) = [];
if exist('stopAtArg', 'var')
    lastValid = find(ARG<=stopAtArg, 1, 'last');
    ARG(lastValid+1:end) = [];
    NEFF(lastValid+1:end) = [];
    clear('stopAtArg');
end;

if infomode
    close(localFigure);
    try figure(oldFigure); catch; end;  %#ok<CTCH>
end;

res = mode;
res.NEFF = NEFF;
res.ARG = ARG;