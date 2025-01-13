function res = modeVsV(mode, region)
% Returns the same mode with V-parameter as the argument

if nargin < 1
    error('%s: at least one mode must be given')
elseif nargin < 2
    region = [];
end

if numel(mode) > 1
    res = [];
    for i = 1:numel(mode)
        r = modeVsV(mode(i), region);
        res = [res r]; %#ok<AGROW>
    end
    return
end

res = mode;
materials = mode.fibreSpec.materials;
if strcmpi(partype(mode), 'wvl')
    lambda = mode.par;
    diameter = mode.ARG;
    res.DIA = mode.ARG;
    res.argtype = 'vpw'; % V-parameter, parameter is wavelength
else % check is done inside partype function
    lambda = mode.ARG;
    diameter = mode.par;
    res.WVL = mode.ARG;
    res.argtype = 'vpd'; % V-parameter, parameter is diameter
end

if strcmpi(region, 'core') && numel(materials) == 3 % core mode in 3LS
    diameter = diameter * mode.fibreSpec.coreCladdingRatio;
    mat1 = materials{1};
    mat2 = materials{2};
elseif strcmpi(region, 'cladding') && numel(materials) == 3 % clad mode in 3LS
    mat1 = materials{2};
    mat2 = materials{3};
else
    mat1 = materials{1};
    mat2 = materials{2};
end

v = diameter * pi ./ (lambda / 1000) .* sqrt(refrIndex(mat1, lambda).^2 - ...
    refrIndex(mat2, lambda).^2);

res.ARG = v;
