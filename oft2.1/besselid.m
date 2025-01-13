function [res] = besselid(n, x)
% Derivative of the modified Bessel function of the first kind

res = besseli(n+1,x) + n./x .* besseli(n,x);
