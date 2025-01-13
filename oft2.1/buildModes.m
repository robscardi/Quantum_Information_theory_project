function [MODES] = buildModes(argument, fibreSpec, modeTask, infomode)
% Calculates the modal curves (n_eff vs d or n_eff vs lambda)
%
% argument  - structure containing type, min and max values
% fibreSpec - fibre specification: materials and core-cladding ratio
% modeTask  - type, nu, harmonic (for arg. type 'wvl'), maxmodes
% infomode  - true|false
% skipmodes - array specifying mode numbers which should not be traced

% Copyright: (cc-by) Subwavelength-diameter fibres team @ Uni Bonn, 2010
% http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com

assert(nargin == 4, 'Wrong input argument number\n');

checkTask(modeTask, fibreSpec);

ModeVsLambda = strcmpi(argument.type, 'WVL');
if ModeVsLambda
    par = modeTask.diameter;
else
    par = modeTask.lambda;
end

%% mode storage
if infomode
    oldFigure = get(0,'CurrentFigure');
    localFigure = figure;
    % ttl = ['\lambda: celltostr(MODETYPE)'];  %#ok<NBRAK>
    
    % Setting for the figure
    hold on;
    if ModeVsLambda
        xlabel('Lambda, nm');
        ylabel('n_{eff}');
        parString = 'radius';
    else
        xlabel('Diameter, \mu{}m');
        ylabel('n_{eff}');
        parString = 'lambda';
    end;
    grid on;
end;

% Check that ARGTYPE is valid
assert(strcmpi(argument.type, 'DIA') || strcmpi(argument.type, 'WVL'), ...
    '%s: ERROR: Wrong mode argument type: %s\n', mfilename, upper(argument.type));

%% SEARCH, TRACE AND STORE MODES
% find modes, trace them, plot them and store them to find intersection
% later

MODES = [];
neffmax = -inf;
neffmin = inf;

tic

if iscell(modeTask.type)
    numberOfModeTypes = length(modeTask.type);
else
    numberOfModeTypes = 1;
end

for i = 1:numberOfModeTypes
    for j = 1:length(modeTask.nu)
        if numberOfModeTypes > 1
            modetype = char(modeTask.type{i});
        else
            modetype = char(modeTask.type);
        end
              
        nu = modeTask.nu(j);
        if ~isfield(modeTask, 'region') && numel(fibreSpec.materials) == 2 % in case of 2LS configuration
            modeTask.region = [];
        end;
        currentTask = struct('modeindex', [nu NaN],...
                'region', modeTask.region,...
                'modetype', modetype,...
                'maxmode', modeTask.maxmode);
        [STARTNEFF, STARTARG, modesDensity] = findModes(argument, par, fibreSpec, currentTask, infomode);
        numberofmodes = length(STARTNEFF);
        if infomode
            fprintf('%s: %g %s mode(s) of order %g found (par=%g, harmonic=%g)\n', mfilename, numberofmodes, upper(char(modetype)), nu, par, NaN);
        end;
        
        %% Trace modes
        
        switch lower(argument.type)
            case 'dia'
                parString = 'lambda';
            case 'wvl'
                parString = 'diameter';
        end
        
        density = max(200, modesDensity);
        if strcmpi(modetype, 'monerie'), density = 500; end % TODO: make it universal for 3LS modes?
        if strcmpi(modetype, 'erdogan'), density = 2000; end % TODO: make it universal for 3LS modes?
        if strncmpi('tsao', modetype, 4), density = 500; end % TODO: make it universal for 3LS modes?
        for m = 1:min(numberofmodes, modeTask.maxmode)
            if infomode
                fprintf('%s: Mode %s (%g,%g) (%s = %5.3f)\n', mfilename, upper(modetype), nu, m, parString, par);
            end;
            
            currentTask.modeindex = [nu m];
            
            if ~isfield(modeTask, 'skipmodes')
                modeTask.skipmodes = [];
            end
            if isempty(find(modeTask.skipmodes == m, 1))
                [newmode] = traceMode(STARTARG(m), STARTNEFF(m), argument, par, ...
                    fibreSpec, currentTask, density, infomode);
            else
                newmode = struct(...
                    'NEFF', STARTNEFF(m), ...
                    'ARG', STARTARG(m), ...
                    'argtype', argument.type, ...
                    'harmonic', argument.harmonic, ...
                    'par', par, ...
                    'fibreSpec', fibreSpec, ...
                    'modeindex', currentTask.modeindex, ...
                    'modetype', currentTask.modetype);
                
                if isfield(currentTask, 'region')
                    newmode.foundin = lower(currentTask.region);
                else
                    newmode.foundin = [];
                end
            end
            MODES = [MODES; newmode]; %#ok<AGROW>
            
            if infomode
                %% Plot mode
                neffmax = max(neffmax, max(newmode.NEFF));
                neffmin = min(neffmin, min(newmode.NEFF));
                if ~ModeVsLambda
                    colour = colourVsLambda(par);
                else
                    colour = 'k';
                end;
                figure(localFigure);
                hold on;
                plot(newmode.ARG, newmode.NEFF, 'Color', colour, 'LineWidth', 2);
                axis([argument.min argument.max 0.9*neffmin neffmax*1.1]);
                drawnow;
            end; % plot modes
        end; % mode_count
    end % correctModeOrder
end % modeType

tElapsed = toc;

if infomode
    fprintf('%s: Modes traced: %g\n', mfilename, length(MODES));
    
    fprintf('Time: %g\n', tElapsed);
    if nargin > 0
        close(localFigure);
    end;
end;

try figure(oldFigure); catch; end; %#ok<CTCH>