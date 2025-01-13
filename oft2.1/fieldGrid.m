function [FG] =  fieldGrid(F)
% Two-layer grid uniform in R and PHI, for displayField2

r1 = F.dr(1) / 2: F.dr(1): F.diam(1) / 2;
r2 = F.diam(1) / 2 + F.dr(2) / 2: F.dr(2): F.diam(2) / 2;
    phi1 = 0: F.dphi(1): 2 * pi - F.dphi(1) / 2;
    phi2 = 0: F.dphi(2): 2 * pi - F.dphi(2) / 2;

[R1, PHI1] = meshgrid(r1, phi1);
[R2, PHI2] = meshgrid(r2, phi2);

ds1 = R1 * F.dr(1) * F.dphi(1);
ds2 = R2 * F.dr(2) * F.dphi(2);

FG = struct('R1', R1, 'R2', R2, 'PHI1', PHI1, 'PHI2', PHI2, 'ds1', ds1, 'ds2', ds2, 'dr', F.dr, 'dphi', F.dphi);