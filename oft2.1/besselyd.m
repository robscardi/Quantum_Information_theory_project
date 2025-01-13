function [res] = besselyd(n, x)
% Derivative of the Bessel function of the second kind

res = 0.5 * (bessely(n-1,x) - bessely(n+1,x));
