function [quadrature, phase] = optical_channel_func(lo,beta, span, symbols, maximum_field,symbol_rate, samples_per_symbol, sample_num, eta, dispesion_array, obs_time)

assert(span> 0 && bitand(span, span - 1) == 0,"span is not a power of two");

%% UNIT OF MEASURMENT
um = 1e-6;

%% CONSTANTS
B = symbol_rate;           %Symbol rate
c = 299792458;

h = 6.62607015e-34;


fiber_width = 9*um;

Aeff = pi*(fiber_width/2).^2;
eps0 = 8.8541878188e-12;
eps = (1.46).^2*eps0;

fc = c/lo.lambda; 

phot_energy = h*fc;

%% PULSE SHAPING

T = 1/B;
sps = samples_per_symbol;
ts = T/sps;
base_impulse = rcosdesign(beta, span, sps);

base_impulse_normalized = base_impulse;



%signal = zeros(span+1, length(t));
Nsample = sps*span +1;
total_sign = zeros(1, Nsample);

for j = 0:(span)
    total_sign(j*sps + 1) = symbols(j+1);
    %signal(j+1, j*sps +1 ) = symbols(j+1);
    %signal(j+1, :) = conv(base_impulse, signal(j+1, :), "same");
end

total_sign = conv(base_impulse/max(base_impulse), total_sign, "same");

%energy_begin = trapz(abs(total_sign).^2)

%% APPLY DISPERSION

dispersed_signal = conv(total_sign, dispesion_array, "same");

%dispersed_energy = trapz(abs(dispersed_signal).^2)
%% FINAL FILTERING

filtered = conv(dispersed_signal, base_impulse_normalized, "same");

%energy_final = trapz(abs(filtered).^2)

%% PHASE NOISE
f = (-Nsample/2:Nsample/2) / ts/Nsample;

D = lo.field*10^(lo.PSD/10);
S_f = D^2./ ((1 + (f / lo.linewidth).^2)*lo.linewidth*2*pi);

white_noise_real = randn(1, (Nsample)+1);
white_noise_imag = randn(1, (Nsample)+1);

white_noise = white_noise_real;
white_noise_ff = fftshift(fft(white_noise));
delta_f_freq_domain = white_noise_ff .* (sqrt(S_f) + sqrt(D)*1i*white_noise_imag);
delta_f_time_domain = ifftshift(ifft(fftshift(delta_f_freq_domain)));
delta_f_time_domain = real(delta_f_time_domain);

phi = cumsum(delta_f_time_domain);

%% PHOTODETECTOR
    middle = sps*span/2 +1;
    half_sample_num = ceil(sample_num/2);
    symbol_sample = filtered(middle-half_sample_num: middle+half_sample_num);
    phi = phi(middle-half_sample_num:middle+half_sample_num);
    %% IN PHASE
    lo.laser = lo.field*exp(1i*phi);
    avg_p = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field + lo.laser).^2)*ts/phot_energy/obs_time;
    avg_m = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field - lo.laser).^2)*ts/phot_energy/obs_time;

    k_p = poissrnd(avg_p);
    k_m = poissrnd(avg_m);

    I_phot = k_p - k_m;

    %% IN QUADRATURE

    lo.laser = lo.field*exp(1i*(phi + pi/2));
    avg_p = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field + lo.laser).^2)*ts/phot_energy/obs_time;
    avg_m = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field - lo.laser).^2)*ts/phot_energy/obs_time;

    k_p = poissrnd(avg_p);
    k_m = poissrnd(avg_m);

    Q_phot = k_p - k_m;


[quadrature, phase] = deal(I_phot, Q_phot);
end 