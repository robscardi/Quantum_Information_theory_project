function result = modeCheckValue(modes)
% Check-sum for the modes array.

intFactor = 1e6;
r = 0;
for i = 1:length(modes)
    mode = modes(i);
    for j = 1:length(mode.NEFF)
        n = uint32(mode.NEFF(j) * intFactor); % make uint our of double
        if strcmpi(mode.argtype, 'DIA')
            a = mode.ARG(j);
        else
            a = mode.ARG(j);
        end;
        a = uint32(a * intFactor); % make uint our of double
        r = bitxor(r, n);
        r = bitxor(r, a);
    end
end;

result = r;
        