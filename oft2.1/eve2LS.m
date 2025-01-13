function [result, lhs, rhs] = eve2LS(d, neff, lambda, fibreSpec, modeTask)
% Eigen-value equations for vector two-layer structure modes

assert(nargin == 5,'%s: Invalid number of parameters', mfilename);

% Assert mode type for 2-layer structure
switch lower(modeTask.modetype)
    case {'hybrid', 'te', 'tm'}
        % ok
    otherwise
        error('Wrong mode type');
end

nu = modeTask.modeindex(1); 

if numel(fibreSpec.materials) == 3
    switch lower(modeTask.region)
        case 'core'
            outerLayer = 2;
        case 'cladding'
            outerLayer = 3;
        otherwise
            error('Wrong region specification');
    end
elseif numel(fibreSpec.materials) == 2
    outerLayer = 2;
else
    error('Wrong number of materials');
end
 
n_inn = refrIndex(fibreSpec.materials{outerLayer - 1}, lambda);
n_out = refrIndex(fibreSpec.materials{outerLayer}, lambda);

if outerLayer == 3
    % Calculate cladding mode of three-layer structure
else 
    if numel(fibreSpec.materials) == 3
        % Convert from cladding diameter to core
        d = d * fibreSpec.coreCladdingRatio;
    end
    % If only two materials are specified, then there is only one diameter
end
    
r = d / 2;

lambda = lambda / 1000; % convert to microns to keep values close to 1 

k0 = 2*pi ./ lambda;
beta = neff .* k0;

betamax = k0 .* n_inn; % max allowed beta value
betamin = k0 .* n_out; % min possible beta value
betamiddle = k0 .* sqrt(n_inn.^2 - n_out.^2); % auxiliary constant for functions

% auxiliary functions from [1, p. 1027]
U = r .* sqrt(betamax.^2 - beta.^2);
W = r .* sqrt(beta.^2 - betamin.^2);
V = r * betamiddle;

% Governing equation F(r, beta, nu) = 0

modetype = char(modeTask.modetype); % conversion from possible cell to string
modetype = modetype(1,:);

switch lower(modetype)
    case 'hybrid' % HE and EH modes
        if nu < 1
            ME = MException('FibreMode:InvalidHYBModeOrder', 'Mode order %g < 1 and mode type is HYB\n', nu);
            throw(ME);
        else
            lhs = (besseljd(nu, U) ./ (U .* besselj(nu, U)) +...
                besselkd(nu, W) ./ (W .* besselk(nu, W))) .*...
                (besseljd(nu, U) ./ (U .* besselj(nu, U)) +...
                n_out.^2 .* besselkd(nu, W) ./...
                (n_inn.^2 .* W .* besselk(nu, W)));
            rhs = (nu * beta ./ (k0 .* n_inn)).^2 .* (V ./ (U .* W)).^4;
            result = lhs - rhs;
                
        end;
    case 'te' %TE modes (Eq. (4)) (added by DP on 2009-01-07)
        if nu ~= 0
            ME = MException('FibreMode:InvalidTEModeOrder', 'Mode order %g ~= 0 and mode type is TE\n', nu);
            throw(ME);
        else
            lhs = besselj(1, U) ./ (U .* besselj(0, U));
            rhs = -besselk(1, W) ./ (W .* besselk(0, W));
            result = lhs - rhs;
        end;
    case 'tm'  % TM mode (Eq. (5)) (added by DP on 2009-01-08)
        if nu ~= 0
            ME = MException('FibreMode:InvalidTMModeOrder', 'Mode order %g ~= 0 and mode type is TM\n', nu);
            throw(ME);
        else
            lhs = n_inn.^2 .* besselj(1, U) ./ (U .* besselj(0, U));
            rhs = -n_out.^2 .* besselk(1, W) ./(W .* besselk(0, W));
            result =  lhs - rhs;
        end;
    otherwise
        ME = MException('FibreMode:InvalidModeType', 'Unknown mode type: %s', upper(modetype));
        throw(ME);
end;