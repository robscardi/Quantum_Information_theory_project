function [result] = normalizePower(F, P)
% Returns field structure normalised to the given power

if nargin == 1 % desired power NOT provided
    P = 1;
end;
factor = fieldPower(F) / P;
F.E1 = F.E1 / sqrt(factor);
F.E2 = F.E2 / sqrt(factor);
F.H1 = F.H1 / sqrt(factor);
F.H2 = F.H2 / sqrt(factor);

result = F;