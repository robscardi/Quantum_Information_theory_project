function [pmpoint] = phaseMatch(mode1, mode2, infomode)
% Tries to find a phase-match point on a between the two given modes
% Returns a PMP structure with quality parameter defines like this:
%   0  - no intersection found
%   1  - fine intersection found
%   -1 - rough intersection found
%
% EXAMPLE
%   See aaGrubskySavchenko.m for example
% 
%
% See also struct_pmp, struct_mode, aagrubskysavchenko, contents

% Copyright: see contents

% NOTE: This function uses the division by two algorithm in the ROUGH
% SEARCH. It seems that it would be better to put new search point to x0
% rather than to (max+min)/2. But this is not correct! Though faster, this
% method may fail in some cases. May be it is possible to make it survive,
% but for now I don't see how to do it and the current method is reliable
% and very fast anyway. Further improvement of rough search will not bring
% improvement in overall performance, since fine search is the most time
% consuming function. KK 2009-06-28

%% Initialization, input check, settings
if nargin < 2
    throw(MException('PhaseMatch:NotEnoughInputArguments', '%s: Not enough input arguments\n', upper(mfilename)));
end;

if nargin < 3
    infomode = true;
end;

argtype = mode1.argtype;
assert(strcmpi(argtype, 'DIA') || strcmpi(argtype, 'WVL'), ...
    '%s: ERROR: Invalid mode argument type: %s\n', mfilename, upper(argtype));
ModeVsLambda = strcmpi(argtype, 'WVL');
assert(ModeVsLambda == strcmpi(mode2.argtype, 'WVL'), ...
    'Both modes should have equal argument type');

assert(ModeVsLambda || (mode1.harmonic == 1 && mode2.harmonic == 1));

grey = [0.5, 0.5, 0.5];

if infomode
    % create a figure for monitoring
    oldFigure = get(0,'CurrentFigure');
    localFigure = figure; hold on;
    
    mode1handle = plot(mode1.ARG, mode1.NEFF, '-o', 'Color', grey, 'MarkerSize', 2, 'MarkerFaceColor', grey);
    mode2handle = plot(mode2.ARG, mode2.NEFF, '-o', 'Color', grey, 'MarkerSize', 2, 'MarkerFaceColor', grey);
    drawnow;
end;

% if strcmpi(argtype, 'WVL')
% The argument is wavelength, so the modal curves can intersect either
% if they correspond to different radii or modes or if they are two
% different harmonics.
%     if par1 ~= par2
%         harmonic = 1;
%     end;

% Leftmost point where both modes are defined:
arg_min = max(mode1.ARG(1), mode2.ARG(1));
% Rightmost point where both modes are defined:
arg_max = min(max(mode1.ARG), max(mode2.ARG));
% Numbers of first(A and C) and last(B and D) points to consider:
% first point of each MODE, where iA,iC >=arg_min
% last point of each MODE, where iB,iD <= r_max
iA = find(mode1.ARG > arg_min, 1, 'first') - 1;
% iB = find(mode1.ARG >= arg_max, 1, 'last');
iB = find(mode1.ARG >= arg_max, 1, 'first');
iC = find(mode2.ARG > arg_min, 1, 'first') - 1;
iD = find(mode2.ARG >= arg_max, 1, 'first');

% Segment of the first mode
A = [mode1.ARG(iA) mode1.NEFF(iA)];
B = [mode1.ARG(iB) mode1.NEFF(iB)];
% Segment of the second mode
C = [mode2.ARG(iC) mode2.NEFF(iC)];
D = [mode2.ARG(iD) mode2.NEFF(iD)];

if infomode
    set(0,'CurrentFigure',localFigure);
    pointAHandle = plot(A(1), A(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
    pointBHandle = plot(B(1), B(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
    pointCHandle = plot(C(1), C(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    pointDHandle = plot(D(1), D(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
end;

%% STEP 1: check for existence of an intersection and harmonic check
if ~IntersectionExists
    if infomode, fprintf('%s: no intersection.\n', mfilename); end;
    exitcode = 0;
    neff1 = NaN;
    neff2 = NaN;
    arg = NaN;
    pmpoint = returnpmpoint;
else % intersection exists
    if infomode, fprintf('%s: intersection exists... ', mfilename); end;
    exitcode = -2;
    searchCount = 0;
    %% STEP 2: Rough search for intersection
    while exitcode < 0
        try
            [iA, iC, x0, y0] = RoughIntersection;
            exitcode = -1;
        catch  %#ok<CTCH>
            if infomode, fprintf('%s: fast rough search failed, trying brute force... ', mfilename); end
            try 
                [iA, iC, x0, y0] = BruteForceRoughSearch;
                exitcode = -1; % brute-force search succeeded
            catch %#ok<CTCH> %
                exitcode = 0; % brute-force search failed---no intersection in fact
            end;
        end;
        
        if exitcode == -1 % intersection exists and rough solution found
            if infomode, fprintf('%s: rough found... ', mfilename); end
            exitcode = -1; % rough intersection found
            neff1 = y0;
            neff2 = y0;
            arg = x0;
            pmpoint = returnpmpoint;
        
            %% STEP 3: Fine search for match point
            [x0, y0, exitcode] = fineSearch;
            if exitcode == 1
                if infomode, fprintf('%s: fine found.\n', mfilename); end
                neff1 = y0;
                neff2 = y0;
                arg = x0;
                pmpoint = returnpmpoint;
            else
                if infomode, fprintf('%s: fine failed, trying rough again... ', mfilename); end
            end;
            searchCount = searchCount + 1; 
            if searchCount == 3
                fprintf('%s: no fine intersection could be found.\n', upper(mfilename)); 
                exitcode = -1;
                pmpoint = returnpmpoint;
                return
            end;
        else % there is no intersection
            if infomode, fprintf('%s: no intersection after brute-force rough search.\n', mfilename); end;
            exitcode = 0;
            neff1 = NaN;
            neff2 = NaN;
            arg = NaN;
            pmpoint = returnpmpoint;
        end;
    end;
end;

try close(localFigure); catch; end %#ok<CTCH>
try set(0,'CurrentFigure',oldFigure); catch; end; %#ok<CTCH>

%% Function: intersection exists
    function functionResult = IntersectionExists
        if infomode
            set(0,'CurrentFigure',localFigure);
            xlim([min(A(1), C(1)) max(B(1), D(1))]);
            ylim([min([A(2), B(2), C(2), D(2)]) max([A(2), B(2), C(2), D(2)])]);
            segment1handle = plot([A(1) B(1)]*2, [A(2) B(2)], 'r');
            segment2handle = plot([C(1) D(1)]*2, [C(2) D(2)], 'r');
            drawnow;
        end;
        
        [segIntersExitCode, x0, y0] = SegmentsIntersection(A, B, C, D, false);
        
        if segIntersExitCode < 1 % no intersection exists
            functionResult = false;
        elseif segIntersExitCode == 1 % intersection of the mode curves exists
            if infomode % show segments intersection point
                set(0,'CurrentFigure',localFigure);
                intersectionPointHandle = plot(x0, y0, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
            end;
            functionResult = true;
            try delete(segment1handle, segment2handle, intersectionPointHandle); catch; end;
        end
    end % function IntersectionExits

%% Function: rough search
    function [iAToReturn, iCToReturn, x0, y0] = RoughIntersection
        roughCounter = 0;
        while iB - iA > 1 || iD - iC > 1 % converge to single segments in both modes
            if infomode
                argmaxLineHandle = plot([arg_max arg_max], [0 2], 'r--');
                argminLineHandle = plot([arg_min arg_min], [0 2], 'r--');
            end;
            arg_mid = (arg_max + arg_min)/2; % middle point
            if infomode
                set(0,'CurrentFigure',localFigure);
                argmidLineHandle = plot([arg_mid arg_mid], [0 2], 'r--');
            end;
            
            if iB - iA > 1 % there is a point between A and B in mode1
                iE = find(mode1.ARG >= arg_mid, 1, 'first') - 1;
            else % A and B are neighbours, take B as E
                iE = iB;
            end;
            %assert(iE > iA);
            E = [mode1.ARG(iE) mode1.NEFF(iE)];
            if infomode
                pointEHandle = plot(E(1), E(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'w');
            end;
            
            if iD - iC > 1 % there is a point between C and D in mode2
                iF = find(mode2.ARG >= arg_mid, 1, 'first') - 1;
            else % C and D are neighbours, take D as F
                iF = iD;
            end;
            %assert(iF > iC);
            F = [mode2.ARG(iF) mode2.NEFF(iF)];
            if infomode
                pointFHandle = plot(F(1), F(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w');
            end;
            
            if infomode
                try delete(argmidLineHandle); catch; end; %#ok<CTCH>
                segmentAEhandle = plot([A(1) E(1)], [A(2) E(2)], 'g', 'LineWidth', 2);
                segmentCFhandle = plot([C(1) F(1)], [C(2) F(2)], 'b', 'LineWidth', 2);
            end;
            [segIntersExitCode, x0, y0] = SegmentsIntersection(A, E, C, F, false);
            
            try delete(argminLineHandle, argmaxLineHandle); catch; end;
            if segIntersExitCode == 1 % AE intersects CF: intersection to the left from the middle point
                %                 arg_max = arg_mid; % move max to the middle
                if iA < iE % only if there is still space to decrease distance
                    iB = iE;
                    B = [mode1.ARG(iB) mode1.NEFF(iB)];
                end;
                if infomode
                    try delete(pointBHandle); catch; end;
                    pointBHandle = plot(B(1), B(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
                end;
                if iC < iF % only if there is still space to decrease distance
                    iD = iF;
                    D = [mode2.ARG(iD) mode2.NEFF(iD)];
                end;
                if infomode
                    try delete(pointDHandle); catch; end;
                    pointDHandle = plot(D(1), D(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
                end;
                arg_max = min(B(1), D(1));
                if infomode
                    try delete(segmentAEhandle, segmentCFhandle); catch; end; %#ok<CTCH>
                end;
            else % AE doesn't intersect CF: intersection to the right from the middle point
                try 
                    assert(SegmentsIntersection(E, B, F, D, false) == 1); % check it; if not true, fast rough search will fail 
                catch ME
                    try delete(segmentAEhandle, segmentCFhandle); catch; end; %#ok<CTCH>
                    rethrow ME
                end;
                if iB > iE % only if there is still space to decrease distance
                    iA = iE;
                    A = [mode1.ARG(iA) mode1.NEFF(iA)];
                end;
                if infomode
                    try delete(pointAHandle); catch; end;
                    pointAHandle = plot(A(1), A(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
                end;
                
                if iD > iF % only if there is still space to decrease distance
                    iC = iF;
                    C = [mode2.ARG(iC) mode2.NEFF(iC)];
                end;
                if infomode
                    try delete(pointCHandle); catch; end;
                    pointCHandle = plot(C(1), C(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
                end;
                arg_min = max(A(1), C(1));
                if infomode
                    try delete(segmentAEhandle, segmentCFhandle); catch; end; %#ok<CTCH>
                end;
            end;
            try delete(pointEHandle, pointFHandle); catch; end; %#ok<CTCH>
            if infomode
                xlim([min(A(1), C(1)) max(B(1), D(1))]);
                ylim([min([A(2), B(2), C(2), D(2)]) max([A(2), B(2), C(2), D(2)])]);
                drawnow;
            end;
            roughCounter = roughCounter + 1; assert(roughCounter < 20);
        end; % while
        iAToReturn = iA;
        iCToReturn = iC;
        % NOTE: The returned intersection point is between AE and CF, even
        % if the actual intersecting segments were EB and FD. The error is
        % negligible, since E is almost on AB and F is almost on CD. The fine
        % search is in any case supposed to find the corrected solution. IF
        % YOU WANT TO REMOVE EVEN THIS ERROR, calculate intersection
        % between AB and CD again.
    end % function RoughIntersection

%% Function: brute-force rough search
    function [iAToReturn, iCToReturn, x0, y0] = BruteForceRoughSearch
        iA = 1;
        iC = 1;
        
        segInters = -1;
        while segInters ~= 1
            A = [mode1.ARG(iA) mode1.NEFF(iA)];
            B = [mode1.ARG(iA+1) mode1.NEFF(iA+1)];
            C = [mode2.ARG(iC) mode2.NEFF(iC)];
            D = [mode2.ARG(iC+1) mode2.NEFF(iC+1)];
            if infomode
                set(0,'CurrentFigure',localFigure);
                try delete(pointAHandle,pointBHandle,pointCHandle,pointDHandle); catch; end; %#ok<CTCH>
                pointAHandle = plot(A(1), A(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
                pointBHandle = plot(B(1), B(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
                pointCHandle = plot(C(1), C(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
                pointDHandle = plot(D(1), D(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
                xlim([min(A(1), C(1)) max(B(1), D(1))]);
                ylim([min([A(2), B(2), C(2), D(2)]) max([A(2), B(2), C(2), D(2)])]);
                drawnow;
            end;
            
            [segInters, x0, y0] = SegmentsIntersection(A, B, C, D, false);
            if segInters ~= 1
                if (B(2) < max(C(2), D(2)) || B(2) > min(C(2), D(2))) && B(1) < D(1)
                    iA = iA + 1;
                else
                    iC = iC + 1;
                end;
            end;
            assert(iA+1 <= length(mode1.ARG) && iC+1 <= length(mode2.ARG));
        end;
        if ~segInters
            assert(1==0);
        end;
        iAToReturn = iA;
        iCToReturn = iC;
    end

%% Function: fine search
% Mode ?zero-functions?--functions which are zero along the MODE curve

% Search for exact intersection point provided the intersecting mode curves
% segments from the previous rough search

    function [x0, y0, exitcodeToReturn] = fineSearch % uses iA and iC from prev. calculations
        convergenceFactor = 1e-7;
%         DeltaFactor = 0.3; % the bigger, the further the fzero margins are off-set from the assumed mode
        
        A = [mode1.ARG(iA) mode1.NEFF(iA)];
        B = [mode1.ARG(iA+1) mode1.NEFF(iA+1)];
        C = [mode2.ARG(iC) mode2.NEFF(iC)];
        D = [mode2.ARG(iC+1) mode2.NEFF(iC+1)];
        
        if ModeVsLambda
            F1 = @(neff, lambda) FibreMode(mode1.par, neff, mode1.modeindex(1), lambda/mode1.harmonic, mode1.materials, mode1.modetype);
            F2 = @(neff, lambda) FibreMode(mode2.par, neff, mode2.modeindex(1), lambda/mode2.harmonic, mode2.materials, mode2.modetype);
        else
            F1 = @(neff, rad) FibreMode(rad, neff, mode1.modeindex(1), mode1.par, mode1.materials, mode1.modetype);
            F2 = @(neff, rad) FibreMode(rad, neff, mode2.modeindex(1), mode2.par, mode2.materials, mode2.modetype);
        end;
        
        % Initial values to enter the WHILE loop:
        [segIntersExitCode, x0, y0] = SegmentsIntersection(A, B, C, D, false);
        fineCount = 0; 
        converged = false;
        while ~converged && segIntersExitCode == 1
            if infomode
                set(0,'CurrentFigure',localFigure);
                try delete(argLineHandle); catch; end; %#ok<CTCH>
                argLineHandle = plot([x0 x0], [0 2], 'r--');
                xlim([min(A(1), C(1)) max(B(1), D(1))]);
                ylim([min([A(2), B(2), C(2), D(2)]) max([A(2), B(2), C(2), D(2)])]);
            end;
            % Find exact betas (for two MODEs) at the found approximate
            % intersection point
            g1 = @(b) F1(b, x0);
            g2 = @(b) F2(b, x0);
            
            % FIRST MODE
            try
                %                     if (sign(g1(A(2))) == sign(g1(B(2))))
                %                         delta = ((B(2)) - (A(2))) * DeltaFactor;
                %                         %phaseMatchFzeroCatchCounter = phaseMatchFzeroCatchCounter + 1;
                %                         if infomode
                %                             x = linspace(C(2) - delta, D(2) + delta, 100);
                %                             OldFigureHandle = gcf;
                %                             figure;
                %                             hold on;
                %                             plot(x,g1(x),'.');
                %                             plot(x,g1(x),'-');
                %                             plot([A(2) A(2)], [g1(A(2)) g1(B(2))], 'r', 'LineWidth', 2);
                %                             plot([B(2) B(2)], [g1(A(2)) g1(B(2))], 'r', 'LineWidth', 2);
                %                             grid on;
                %                             figure(OldFigureHandle);
                %                         end;
                %                     else
                %                         delta = 0;
                %                     end;
                %                     %                 neff1 = fzero(g1, [A(2) - delta, B(2) + delta], options);
                neff1 = findRoot(g1, [A(2), B(2)], (x0-min(A(1), B(1)))/(abs(A(1)-B(1))));
                E = [x0, neff1];
                newARG = [mode1.ARG(1:iA) E(1) mode1.ARG(iA+1:end)];
                newNEFF = [mode1.NEFF(1:iA) E(2) mode1.NEFF(iA+1:end)];
                mode1.ARG = newARG; mode1.NEFF = newNEFF;
                if infomode
                    delete(mode1handle);
                    mode1handle = plot(mode1.ARG, mode1.NEFF, '-o', 'Color', grey, 'MarkerSize', 2, 'MarkerFaceColor', grey);
                    try delete(pointEHandle); catch; end; %#ok<CTCH>
                    pointEHandle = plot(E(1), E(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'w');
                end;
            catch ME %#ok<NASGU>
                break;
            end;
            
            % SECOND MODE
            try
                neff2 = FindRoot(g2, [C(2), D(2)], (x0-min(C(1), D(1)))/(abs(D(1)-C(1))));
                F = [x0, neff2];
                newARG = [mode2.ARG(1:iC) F(1) mode2.ARG(iC+1:end)];
                newNEFF = [mode2.NEFF(1:iC) F(2) mode2.NEFF(iC+1:end)];
                mode2.ARG = newARG; mode2.NEFF = newNEFF;
                if infomode
                    delete(mode2handle);
                    mode2handle = plot(mode2.ARG, mode2.NEFF, '-o', 'Color', grey, 'MarkerSize', 2, 'MarkerFaceColor', grey);
                    try delete(pointFHandle); catch; end; %#ok<CTCH>
                    pointFHandle = plot(F(1), F(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w');
                end;
            catch ME %#ok<NASGU>
                break;
            end;
            
            [segIntersExitCode, x0, y0] = SegmentsIntersection(A, E, C, F, false);
            if segIntersExitCode == 1 % intersection seems to be to the left from E and F
                B = E;
                D = F;
            else % intersection is not to the left from E and F
                [segIntersExitCode, x0, y0] = SegmentsIntersection(E, B, F, D, false);
                if segIntersExitCode == 1 % intersection seems to be to the right from E and F
                    iA = iA + 1;
                    iC = iC + 1;
                    A = E;
                    C = F;
                end; % intersection lost; repeat rough search with new points in the modes
            end;
            if infomode
                try delete(pointAHandle); catch; end; %#ok<CTCH>
                pointAHandle = plot(A(1), A(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
                try delete(pointCHandle); catch; end; %#ok<CTCH>
                pointCHandle = plot(C(1), C(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
                try delete(pointBHandle); catch; end; %#ok<CTCH>
                pointBHandle = plot(B(1), B(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
                try delete(pointDHandle); catch; end; %#ok<CTCH>
                pointDHandle = plot(D(1), D(2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            end;
            %                 try
            %                     if sign(g2(C(2))) == sign(g2(D(2)))
            %                         delta = ((D(2)) - (C(2))) * DeltaFactor;
            %                         %phaseMatchFzeroCatchCounter = phaseMatchFzeroCatchCounter + 1;
            %                         if infomode
            %                             x = linspace(C(2) - delta, D(2) + delta, 100);
            %                             OldFigureHandle = gcf;
            %                             figure;
            %                             hold on;
            %                             plot(x,g1(x),'.');
            %                             plot(x,g1(x),'-');
            %                             plot([C(2) C(2)], [g1(C(2)) g1(D(2))], 'r', 'LineWidth', 2);
            %                             plot([D(2) D(2)], [g1(C(2)) g1(D(2))], 'r', 'LineWidth', 2);
            %                             grid on;
            %                             figure(OldFigureHandle);
            %                         end;
            %                     else
            %                         delta = 0;
            %                     end;
            %
            %                     %                 neff2 = fzero(g2, [C(2) - delta, D(2) + delta], options);
            %                     neff2 = FindRoot(g2, [C(2) - delta, D(2) + delta]);
            %
            %                 catch ME %#ok<NASGU>
            %                     % If there is no zero point between C(2) and D(2), one of
            %                     % this points is to far away from the real zero point. In
            %                     % other words, the tolerance of fzero is to smal. To get
            %                     % rid of this error increase DeltsFactor.
            %                     %                 rethrow(ME);
            %                     neff1 = y0;
            %                     neff2 = y0;
            %                     arg1 = x0;
            %                     arg2 = x0;
            %                     assert(exitcode == -1);
            %                     pmpoint = returnpmpoint;
            %                     return; % returns rough result
            %                 end;
            
            %                 BOLD = B; B = [arg, neff1];
            %                 DOLD = D; D = [arg, neff2];
            %                 AOLD = A; %#ok<NASGU>
            %                 COLD = C; %#ok<NASGU>
            
            %             elseif segIntersExitCode == 2
            %                 fprintf('%s: ERROR: segments coincide\n', mfilename);
            %
            %             elseif segIntersExitCode == -1 % intersection lies only on segment AB, so shift CD
            %                 if x0 < C(1)
            %                     D = C;
            %                     iC = iC - 1;
            %                     C(1) = mode2.ARG(iC); C(2) = mode2.NEFF(iC);
            %                 elseif x0 > D(1)
            %                     C = D;
            %                     iD = iD + 1;
            %                     D(1) = mode2.ARG(iD); D(2) = mode2.NEFF(iD);
            %                 else
            %                     throw(MException('PhaseMatch:IntersectionDoesNotMatchExitcode'));
            %                     % segIntersExitCode == -1, but C(1) <= x <= D(1)
            %                 end;
            %
            %             elseif segIntersExitCode == -2 % intersection lies only on segment CD, so shift AB
            %                 if x0 < A(1)
            %                     B = A;
            %                     iA = iA - 1;
            %                     A(1) = mode1.ARG(iA); A(2) = mode1.NEFF(iA);
            %                 elseif x0 > B(1)
            %                     A = B;
            %                     iB = iB + 1;
            %                     B(1) = mode1.ARG(iB); B(2) = mode1.NEFF(iB);
            %                 else
            %                     throw(MException('PhaseMatch:IntersectionDoesNotMatchExitcode'));
            %                     % segIntersExitCode == -2, but A(1) <= x <= B(1)
            %                 end;
            %
            %             elseif segIntersExitCode == -3 % intersection lies neither on segment CD nor on segment AB
            %                 % This condition is needed if the phase match point is not in
            %                 % the section that was checked, but in the other.
            %                 % Here the points are set so the other section can be checked
            %                 A = B; C = D;
            %                 B = BOLD; D = DOLD;
            %
            %             else
            %                 throw(MException('PhaseMatch:InvalidSegmentsIntersectionExitcode'));
            %             end;
            fineCount = fineCount + 1; assert(fineCount < 1000);
            converged = abs(E(2)-F(2)) / (E(2)+F(2)) < convergenceFactor;
        end; % while
        
        if converged
            exitcodeToReturn = 1;
        else
            exitcodeToReturn = -1;
        end;
    end %function fine search

%% Function: returnpmpoint
function res = returnpmpoint
    res = struct('neff1', neff1, 'neff2', neff2, ...
        'arg', arg, 'argtype', argtype, ...
        'mode1', mode1, 'mode2', mode2, ...
        'quality', exitcode);
    % 1 for fine intersection, -1 for rough, -2 for intersection exists
    % but couldn't be found, 0 for no intersection
    try
        res = setHarmonic(res);
    catch ME
        if strcmpi(ME.message, 'Cannot determine the harmonic')
            res.harmonic = NaN;
        else
            rethrow ME;
        end
    end
end

end % main function


function [exitcode, x, y] = SegmentsIntersection(A, B, C, D, infomode)
% exitcode:
%   0 if no interesection point 
%        Lines parallel or one of the segments is a point not lieing on the
%        second line. In this case (x,y) = (NaN,NaN)
%   -1 if intersection point is outside of segment AB
%   -2 if intersection point is outside of segment CD
%   -3 if intersection point is outside of both segments
%   1 if intersection point is found
%   2 if segments lie on the same line (multiple intersection points)
%
% (by) Konstantin Karapetyan, Dimitri Pritzkau 2008/2009
% kotya.karapetyan@gmail.com

exitcode = 0;
x = NaN;
y = NaN;

if nargin < 4
    disp([mfilename 'ERROR: Not enough parameters']);
    return;
elseif nargin == 4
    infomode = true;
end;

if infomode
    oldFigure = get(0,'CurrentFigure');
    localFigure = figure; hold on;
end;

Ax = A(1);
Bx = B(1);
Ay = A(2);
By = B(2);
if infomode
    set(0,'CurrentFigure',localFigure);
    plot([Ax Bx], [Ay By], '-m');
    plot([Ax Bx], [Ay By], '.m');
end;

Cx = C(1);
Dx = D(1);
Cy = C(2);
Dy = D(2);
if infomode
    figure(localFigure);
    plot([Cx Dx], [Cy Dy], 'c');
    plot([Cx Dx], [Cy Dy], '.c'); % points
end;

CDy = (Dy - Cy);
CDx = (Dx - Cx);
ABy = (By - Ay);
ABx = (Bx - Ax);

%% Kotyas version

% numerator = (Ay - Cy) * ABx * CDx - Ax * ABy * CDx + Cx * CDy * ABx;
% denominator = CDy * ABx - ABy * CDx;
% 
% if abs(denominator) < 1E-15
%     x = NaN;
%     y = NaN;
%     if abs(numerator) < 1E-15 % lines coincide
%         exitcode = 2;
%         if infomode
%             disp('SEGMENTSINTERSECTION: Lines coincide -- multiple intersection points');
%         end;
%     else
%         exitcode = 0; % lines are parallel
%         if infomode
%             disp('SEGMENTSINTERSECTION: Segments are parallel');
%         end;
%     end;
%     return;
% end;
% 
% x = numerator / denominator;
% y = Ay + ABy / ABx * (x - Ax);
% 


%% Dimitri's version

% get rid of infinite gradients
if ABx == 0 && CDx == 0 % both lines are vertical
    x = NaN;
    y = NaN;
elseif ABx == 0 % only AB is vertical
    if ABy == 0 % A and B coincide
        if (Ay - Cy) / (Ax - Cx) == CDy / CDx % A & B lie on CD
            x = Ax; y = Ay;
        else
            x = NaN; y = NaN;
        end;
    else 
        m2 = CDy / CDx;
        x = Ax; % Since AB is vertical, take its x...
        y = m2 * (Ax - Cx) + Cy; % ...and find y from CD
    end;
elseif  CDx == 0 % only CD is vertical
    if  CDy == 0 % C and D coincide
        if (Cy - Ay) / (Cx - Ax) == (ABy / ABx)
            x = Cx; y = Cy;
        else
            x = NaN; y = NaN;
        end;
    else
        m1 = ABy / ABx;
        x = Cx;
        y = m1 * (Cx - Ax) + Ay;
    end
else
    m1 = ABy / ABx;
    m2 = CDy / CDx;
    x = Ax + ((m2 * (Ax - Cx) + Cy - Ay ) / (m1 - m2));
    y0 = m1 * (-Ax) + Ay;
    y = m1 * x +y0;
        
    if isinf(x) % both lines are horizontal
        x = NaN; y = NaN;
    end;
end;
if infomode
    set(0,'CurrentFigure',localFigure);
    plot(x, y, '.k', 'MarkerSize', 8);
end;

%% EXITCODE section
if isnan(x)
    if Ax == Cx
        exitcode = 2; % lines coincide
        if infomode
            disp([mfilename ': Lines coincide']);
        end;
    else
        exitcode = 0;
        if infomode
            disp([mfilename ': No intersection: lines are parallel or one segment is an external point for the second']);
        end;
    end;
elseif (x >= min(Ax, Bx)) && (x <= max(Ax, Bx)) && (x >= min(Cx, Dx)) && (x <= max(Cx, Dx))
    exitcode = 1; % the intersection point lies within the segments
    if infomode
        disp([mfilename ': Intersection point found']);
    end;
elseif (x > min(Ax, Bx)) && (x < max(Ax, Bx))
    exitcode = -1; % intersection point is inside the segment AB
    if infomode
        disp([mfilename ': Intersection point lies inside segment AB only']);
    end;
elseif (x > min(Cx, Dx)) && (x < max(Cx, Dx))
    exitcode = -2; % intersection point is inside the segment CD
    if infomode
        disp([mfilename ': Intersection point lies inside segment CD only']);
    end;
else
    exitcode = -3; % the intersection point lies outside both segments
    if infomode
        disp([mfilename ': Intersection point lies outside both segments']);
    end;
end;


% % ONLY FOR TESTING
% % Chek if the intersection point is completely outside the range
% if isnan(x) || isnan(x) ||...
%         x > max([Ax Bx]) || x > max([Cx Dx]) ||...
%         x < min([Ax Bx]) || x < min([Cx Dx]) ||...
%         y > max([Ay By]) || y > max([Cy Dy]) ||...
%         y < min([Ay By]) || y < min([Cy Dy])
%     exitcode = 0; % no intersection
% else
%     exitcode = 1;
% end;
% 
if infomode
    figure(localFigure);
%     plot(Bx, By, '.r'); % points
%     plot(Dx, Dy, '.b'); % points
    grid on;
    fprintf('%s:x = %g, y = %g, exitcode = %g\n', mfilename, x , y, exitcode);
    try close(localFigure); catch; end %#ok<CTCH>
    try figure(oldFigure); catch; end; %#ok<CTCH>
end;

end % segmentsIntersection function





















