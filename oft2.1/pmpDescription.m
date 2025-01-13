function result = pmpDescription(pmp, mode)
% Creates a text description of a phase-matching point

if nargin == 1
    mode = 'long';
end

pmp = setHarmonic(pmp);

if numel(pmp) > 1
    result = cell(1,numel(pmp));
    for i=1:numel(pmp)
        result{i} = pmpDescription(pmp(i), mode);
    end
    return
end

if strcmpi(mode, 'long') 
    result = [sprintf('  %s intersects \n', modeDescription(pmp.mode1))...
        sprintf('  %s\n', modeDescription(pmp.mode2))];
    
    if strcmpi(pmp.argtype, 'dia')
        result = [result sprintf('  at d = %g um\n', pmp.arg)];
    elseif strcmpi(pmp.argtype, 'wvl')
        result = [result sprintf('  at lambda = %g nm\n', pmp.arg)];
    else
        error ('Invalid argtype');
    end;
elseif strcmpi(mode, 'short') 
    result = sprintf('%s to %s', modeDescription(pmp.mode1, false), modeDescription(pmp.mode2, false));
else
    error('INVALID MODE');
end