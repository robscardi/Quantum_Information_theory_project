function [res] = besselkd(n, x)
%  Derivative of the modified Bessel function of the second kind

res = -0.5 * (besselk(n-1,x) + besselk(n+1,x));

