close all
clear variables

t = (-4:0.02:4);
T = 1;

b = rcosdesign(0.5, 8, (length(t)-1)/8);
energy = trapz(b.^2);
figure
plot(t, b)
bb = fft(b);
fs = 1/(0.2);
L = length(t);
f = fs/L*(-L/2:L/2-1);
figure
hold on
plot(f, abs(fftshift(bb)));
plot(f, angle(fftshift(bb)));

h = ifft(bb/max(abs(bb)));
energy = trapz(h.^2);

figure
plot(t, h);

hh = bb.*bb;
h = ifftshift(ifft(hh));
figure
plot(t, h)

b = b*(1*exp(i*pi/4));
bb = fft(b);
fs = 1/(0.2);
L = length(t);
f = fs/L*(-L/2:L/2-1);
figure
hold on
plot(f, abs(fftshift(bb)));
plot(f, angle(fftshift(bb)));

