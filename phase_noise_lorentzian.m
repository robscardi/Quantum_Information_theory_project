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
sample_num = 100;

T = 1/B;
sps = 100;
ts = T/sps;

lo.linewidth = 1*kHz;
lo.PSD = 80;
lo.lambda = 1550*nm;
lo.field = 1e3;
fc = c/lo.lambda;

f = (-sample_num/2:sample_num/2) /ts/sample_num;

D = lo.field;
S_f = D^2./ ((1 + (f / lo.linewidth).^2)*lo.linewidth*2*pi);

white_noise_real = randn(1, sample_num+1);
white_noise_imag = randn(1, sample_num+1);

white_noise = white_noise_real + 1j * white_noise_imag;
delta_f_freq_domain = white_noise .* sqrt(S_f);
delta_f_time_domain = ifft(delta_f_freq_domain);
delta_f_time_domain = real(delta_f_time_domain);

phi = cumsum(delta_f_time_domain);

figure
grid on
stem(phi)

ff = fftshift(fft(phi));
figure
grid on
plot(10*log10(abs(ff)))
