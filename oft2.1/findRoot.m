function result = findRoot(g, LIMITS, p)
% Works on top of fzero, tries to make it work in larger number of cases.
%
% FINDROOT tries to find the root (change-of-sign point) of a function
% g, assuming that this should happen between LIMITS(1) and LIMITS(2). The
% function uses MATLAB's FZERO but can handle some cases when LIMTS are not
% directly applicable to FZERO (the function is not real or there is no
% change of sign in between). 
% 
% The search starts in the center between the LIMITS. This can be changed
% by providing the relative distance from the first limit, P: P == 0 means
% start search from LIMITS(1), P == 0.5 means start search from the center,
% P == 1 means start search from LIMITS(2).
%
% See source code for algorithm.
%
% ALGORITHM
% TODO: describe algorithm

%% Settings
ToleranceFactor = 1e-6;
AttemptMax = 10;
ExtensionFactor = 1.2; % together, on both sides

%% Check input data
assert(length(LIMITS) == 2);
assert(LIMITS(1) ~= LIMITS(2), 'FINDROOT:EqualLimits');

LIMITS = sort(LIMITS);
x1 = LIMITS(1);
x2 = LIMITS(2);
assert(nargin == 2 || nargin == 3);
if nargin == 2
    p = 0.5;
end;

x0 = x1 + (x2 - x1) * p;

%%
attempt = 1;
NofSteps = 5;

while attempt < AttemptMax;
    dx = (x2 - x1) / NofSteps;
    NofStepsExtended = round(NofSteps * ExtensionFactor);
    if mod(NofStepsExtended, 2) == 0 % assurance that NofStepsExtended is odd
        NofStepsExtended = NofStepsExtended + 1;
    end;
    
    ExtensionCount = NofStepsExtended - NofSteps;
    
    X = x1 - (ExtensionCount / 2) * dx: dx: x2 + (ExtensionCount / 2) * dx;
    G = g(X);
    
    StartingInterval = find(X > x0, 1, 'first') - 1;
    for i = 1:NofStepsExtended;
        % Alternating intervals starting from the center: center, center+1,
        % center-1, center+2, center-2 etc.
        interval = StartingInterval + floor(i/2) * (-1)^(i-1);
        if interval < 1
            continue;
        end;
        xLeft = X(interval); xRight = X(interval+1);
        gLeft = G(interval); gRight = G(interval+1);
        
        % Try to find zero
        if isreal([gLeft, gRight]) && sign(gLeft(1)) == -sign(gRight)
            fzerooptions = optimset('TolX', max(10*eps, dx*ToleranceFactor));
            try
                [root, ~, fzerores] = fzero(g, [xLeft, xRight], fzerooptions);
                if fzerores ~= 1
                    continue;
                end;
                result = root;
                return
            catch %#ok<CTCH>
                continue;
                % In principle, we should never get here: interval ends are
                % fine; if the point found is not the root but the 
                % singularity point, it will be intercepted by IF above. So
                % it's better to take care of the case when we are getting 
                % here. (KK 2011-11-02)
            end;
        end;
    end;
    NofSteps = NofSteps * 2 - 1;
    attempt = attempt + 1;
end;

ME = MException('FINDROOT:ReachedMaxNumberOfTrials', ...
    '%s: Reached max value of ATTEMPT (%g) before root was found\n', ...
    upper(mfilename), AttemptMax);
throw(ME);

