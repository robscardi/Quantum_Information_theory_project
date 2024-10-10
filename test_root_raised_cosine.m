close all
clear variables

t = (-4:0.02:4);
T = 1;
%% BASE ROOT RAISED COSINE
base_impulse = rcosdesign(0.5, 8, (length(t)-1)/8);
energy_1 = trapz(base_impulse.^2);
figure
plot(t, base_impulse)
base_impulse_ft = fft(base_impulse);
fs = 1/(0.2);
L = length(t);
f = fs/L*(-L/2:L/2-1);
figure
hold on
plot(f, abs(fftshift(base_impulse_ft)));
plot(f, angle(fftshift(base_impulse_ft)));

%% TEST THE SYMBOL ENERGY AFTER NORMALIZED FOURIER TRANSFORM
base_impulse_normalized = ifft(base_impulse_ft/max(abs(base_impulse_ft)));
energy_2 = trapz(base_impulse_normalized.^2);

figure
plot(t, base_impulse_normalized);

%% TEST SYMBOL ENERGY AFTER SECOND ROOT RAISED COSINE
hh = base_impulse_ft.*base_impulse_ft;
after_filtering = ifftshift(ifft(hh));
figure
plot(t, after_filtering)

% base_impulse = base_impulse*(1*exp(1i*pi/4));
% base_impulse_ft = fft(base_impulse);
% fs = 1/(0.2);
% L = length(t);
% f = fs/L*(-L/2:L/2-1);
% figure
% hold on
% plot(f, abs(fftshift(base_impulse_ft)));
% plot(f, angle(fftshift(base_impulse_ft))); 















