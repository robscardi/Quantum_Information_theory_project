function res = partype(modes)
% Returns the parameter type ('wvl' or 'dia')

if strcmpi(modes(1).argtype, 'dia') || strcmpi(modes(1).argtype, 'vpw')
    res = 'wvl';
elseif strcmpi(modes(1).argtype, 'wvl') || strcmpi(modes(1).argtype, 'vpd')
    res = 'dia';
else
    error('%s: invalid argument type', mfilename);
end

for i = 2:numel(modes)
    if ~strcmpi(partype(modes(i)), partype(modes(1)))
        error('%s: inconsistent parameter types\n', mfilename);
    end
end 