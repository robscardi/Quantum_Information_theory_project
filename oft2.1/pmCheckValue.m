function result = pmCheckValue(pmPoints, withOverlap)

% Copyright: (cc-by) Subwavelength-diameter fibres team @ Uni Bonn, 2010 
% http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com


if nargin == 1 
    withOverlap = true;
end;

intFactor = 1e6;
r = 0;
pmPoints = setHarmonic(pmPoints);
for i = 1:length(pmPoints)
    pm = pmPoints(i);
    n = uint32(pm.neff1 * intFactor); % make uint our of double
    a = uint32(pm.arg * intFactor); % make uint our of double
    if pm.arg > 50 % argument is in nanometers
        a = uint32(pm.arg * intFactor / 1000); % reduce by 1000 to keep numbers for XOR at the same order
    end
    if pm.quality == -1 % rough result, consider only the first digits of it
        n = uint32(n / 1000);
        a = uint32(a / 1000);
    end;
    if withOverlap
        o = uint32(overlapIntegral(pm) * intFactor); 
        % Note: for very small overlaps the value is calculated wrongly. 
        % These values will be just ignored by uint32.
    else
        o = uint32(0); % don't use overlap (for speed)
    end;
    r = bitxor(r, n);
    r = bitxor(r, a);
    r = bitxor(r, o);
end;

result = r;
        