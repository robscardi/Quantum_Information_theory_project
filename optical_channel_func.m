function [quadrature, phase] = optical_channel_func(lo, span, symbol, maximum_field,symbol_rate, samples_per_symbol, eta, dispesion_array, obs_time)

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
beta = 0.5;
sps = samples_per_symbol;
ts = T/sps;
base_impulse = rcosdesign(beta, span, sps);



base_impulse_normalized = base_impulse/(max(abs(fft(base_impulse))));

symbol_inphase = 1*rand(1,span+1);
symbol_inquadrature = 1*rand(1, span+1);
symbol_inphase(span/2 +1 ) = real(symbol);
symbol_inquadrature(span/2 +1) = imag(symbol); 


t = linspace(-((span+1)/2)*T,  ((span+1)/2)*T, sps*span +1);

signal = zeros(span+1, length(t));
total_sign = zeros(1, length(t));

for j = 0:(span)
    total_sign(j*sps + 1) = symbol_inphase(j+1)+1i*symbol_inquadrature(j+1) ;
    signal(j+1, j*sps +1 ) = symbol_inphase(j+1)+1i*symbol_inquadrature(j+1);
    signal(j+1, :) = conv(base_impulse, signal(j+1, :), "same");
end

total_sign = conv(base_impulse_normalized, total_sign, "same");

%energy_begin = trapz(abs(total_sign).^2)

%% APPLY DISPERSION

dispersed_signal = conv(total_sign, dispesion_array, "same");

%dispersed_energy = trapz(abs(dispersed_signal).^2)
%% FINAL FILTERING

filtered = conv(dispersed_signal, base_impulse_normalized, "same");

%energy_final = trapz(abs(filtered).^2)

%% PHASE NOISE

fmin = fc - lo.linewidth/2;
fmax = fc + lo.linewidth/2;

D = 10^(lo.PSD/10);
var = D*((1/fmin)-(1/fmax));
std = sqrt(var);

sample_num = ceil(obs_time/ts);
sample_num = (mod(sample_num, 2) == 0)*(sample_num) + (mod(sample_num,2)==1)*(sample_num+1); 

phi = zeros(1, sample_num+1);
random_walk = std*randn(1, length(phi)-1);

for i = 1:1:length(phi)-1
    phi(i+1) = phi(i) + random_walk(i);
end

%% PHOTODETECTOR
    middle = sps*span/2 +1;
    half_sample_num = ceil(sample_num/2);
    symbol_sample = filtered(middle-half_sample_num: middle+half_sample_num);
    %% IN PHASE
    lo.laser = lo.field*exp(1i*phi);
    avg_p = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field + lo.laser).^2)*ts/phot_energy;
    avg_m = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field - lo.laser).^2)*ts/phot_energy;

    k_p = poissrnd(avg_p);
    k_m = poissrnd(avg_m);

    I_phot = k_p - k_m;

    %% IN QUADRATURE

    lo.laser = lo.field*exp(1i*(phi + pi/2));
    avg_p = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field + lo.laser).^2)*ts/phot_energy;
    avg_m = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field - lo.laser).^2)*ts/phot_energy;

    k_p = poissrnd(avg_p);
    k_m = poissrnd(avg_m);

    Q_phot = k_p - k_m;


[quadrature, phase] = deal(I_phot, Q_phot);
end 