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

B = 400e6;
c = 299792458;

T = 1/B;
sps = 100;
ts = T/sps;
sample_num = sps*16;

fs = 1/ts;

lo.linewidth = 1*kHz;
lo.PSD = -80;
lo.lambda = 1550*nm;
lo.field = 1e1;
fc = c/lo.lambda;

f = (-sample_num/2:sample_num/2) /ts/sample_num;

D = lo.field*10^(lo.PSD/10);
S_f = D^2./ ((1 + (f / lo.linewidth).^2)*lo.linewidth*2*pi);

T_noise = ifft(S_f);

white_noise_real = randn(1, (sample_num)+1);
white_noise_imag = randn(1, (sample_num)+1);

white_noise = white_noise_real;
white_noise_ff = fftshift(fft(white_noise));
delta_f_freq_domain = white_noise_ff .* (sqrt(S_f) + sqrt(D)*1i*white_noise_imag);
delta_f_time_domain = ifftshift(ifft(fftshift(delta_f_freq_domain)));
delta_f_time_domain = real(delta_f_time_domain);


phi = cumsum(delta_f_time_domain);

figure
grid on
plot(abs(white_noise))

figure
grid on
plot(10*log10(abs(S_f)))

figure
grid on
plot(phi, LineWidth=2)
ylabel("Lorentzian random walk", FontSize=40)
xlabel("Sample number", FontSize=40)
grid on
xlim([0 sample_num])
ax = gca; % Get the current axis
ax.FontSize = 30;

ff = fftshift(fft(phi));
figure
plot(f/MHz, 20*log10(abs(ff)/lo.field), LineWidth=2)
ylabel("Phase noise PSD [dBc/Hz]", FontSize=35)
xlabel("Frequency [MHz]", FontSize=35)
grid on
ax = gca; % Get the current axis
ax.FontSize = 30;