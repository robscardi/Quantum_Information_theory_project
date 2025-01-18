close all
clearvars
um = 1e-6;
nm = 1e-9;
pm = 1e-12;
MHz = 1e6;
kHz = 1e3;
km = 1e3;
ps = 1e-12;
mW = 1e-3;

c = 299792458;
sample_num = 10000;

lo.linewidth = 1000000*kHz;
lo.PSD = 40;
lo.lambda = 1550*nm;
lo.field = 1e1;
fc = c/lo.lambda;

fmin = fc - lo.linewidth/2;
fmax = fc + lo.linewidth/2;

D = 10^(lo.PSD/10);
var = D*((1/fmin)-(1/fmax));
std = sqrt(var);

phi = zeros(1, sample_num+1);
random_walk = std*randn(1, length(phi)-1);

for i = 1:1:length(phi)-1
    phi(i+1) = phi(i) + random_walk(i);
end

figure
grid on
plot(phi)

ff = fftshift(fft(phi));
figure
grid on
plot(10*log10(abs(ff/sample_num)))
