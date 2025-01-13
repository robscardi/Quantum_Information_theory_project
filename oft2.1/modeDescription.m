function result = modeDescription(mode, withParameter, withHarmonic)
% Returns text description of the mode(s).
%
% If withParameter is true, the mode parameter (wavelength or diameter) is
% appended in parenthesis.
% If numel(mode)>1, returns a cell array 

if nargin == 1
    withParameter = true;
    withHarmonic = true;
end;

if nargin == 2
    withHarmonic = true;
end;

if numel(mode) > 1
    result = cell(1,numel(mode));
    for i=1:numel(mode)
        result{i} = modeDescription(mode(i), withParameter, withHarmonic);
    end
    return
end

if nargin == 1
    withParameter = true;
end;

if strcmpi(mode.modetype, 'te') || strcmpi(mode.modetype, 'tm')
    % TE or TM + index
    result = [upper(mode.modetype) ' ' mat2str(mode.modeindex(1)) ' ' mat2str(mode.modeindex(2))];
elseif strcmpi(mode.modetype, 'hybrid')
    % HE or EH + index
    if mod(mode.modeindex(2), 2) == 1 % odd modes
        text = 'HE';
        modeNumber = ceil(mode.modeindex(2) / 2);
    else
        text = 'EH';
        modeNumber = floor(mode.modeindex(2) / 2);
    end;

    result = [upper(text) ' ' mat2str(mode.modeindex(1)) ' ' mat2str(modeNumber)];
elseif strncmp(mode.modetype, 'tsao', 4)
    result = ['TSAO' upper(mode.modetype(5:6)) ' ' mat2str(mode.modeindex(1)) ' ' mat2str(mode.modeindex(2))];
elseif strcmpi(mode.modetype, 'erdogan') || strcmpi(mode.modetype, 'monerie')
    result = [upper(mode.modetype) ' ' mat2str(mode.modeindex(1)) ' ' mat2str(mode.modeindex(2))];
else
    error('%s: invalid mode type (%s)\n', upper(mfilename), mode.modetype)
end;

if withHarmonic || withParameter
    result = [result ' ('];
end

% add harmonic
if withHarmonic
    result = [result mat2str(mode.harmonic) ' hrm.'];
end

if withHarmonic && withParameter
    result = [result ', '];
end

% add parameter (diameter or wavelength)
if withParameter
    if strcmpi(partype(mode), 'wvl')
        result = [result mat2str(mode.par, 6) ' nm)'];
    else
        result = [result mat2str(mode.par, 6) ' um)'];
    end;
else
    if withHarmonic
        result = [result ')'];
    end
end