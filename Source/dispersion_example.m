close all
clear variables

um = 1e-6;
nm = 1e-9;
pm = 1e-12;
MHz = 1e6;
kHz = 1e3;
km = 1e3;
ps = 1e-12;
mW = 1e-3;
c = 299792458;
B = 400e6;
span = 4;
n_bit = 3;
M = 2^n_bit;
total_symbols = qammod(0:M-1, M);
symbol_vec = unique(real(total_symbols));
max_symbol = max(symbol_vec);

random_indices = randi(length(total_symbols), 1, span+1);
symbols = total_symbols(random_indices);

beta = 0.2;
T = 1/B;
sps = 100;
ts = T/sps;
base_impulse = rcosdesign(beta, span, sps);

base_impulse_normalized = base_impulse;

Nsample = sps*span +1;
total_sign = zeros(1, Nsample);
signal = zeros(span+1, Nsample);

for j = 0:(span)
    total_sign(j*sps + 1) = symbols(j+1);
    signal(j+1, j*sps +1 ) = symbols(j+1);
    signal(j+1, :) = conv(base_impulse, signal(j+1, :), "same");
end

total_duration = span*T;
fs = 1/ts;
L = span*sps;
f = fs/L*(-L/2:L/2-1);
lambda = 1550*nm;
fc = c/lambda;

lambda_vector = c./(f+fc);
Communication_lenght = 250*km;

D = 17*(ps/(nm*km));
beta = D*((lambda.*f).^2*pi/c);

ff = exp(-1i*beta*Communication_lenght);

tt = ifftshift(ifft(ifftshift(ff)));

total_sign = conv(base_impulse/max(base_impulse), total_sign, "same");
dispersed_signal = conv(total_sign, tt, "same");
filtered = conv(dispersed_signal, base_impulse_normalized, "same");

single_simb_disp = zeros(size(signal));
for i=1:size(signal,1)
    single_simb_disp(i,:) = conv(signal(i,:), tt, "same");
    single_simb_disp(i,:) = conv(single_simb_disp(i,:), base_impulse_normalized, "same");
end

single_simb_undisp = zeros(size(signal));
for i=1:size(signal,1)
    single_simb_undisp(i,:) = signal(i,:);
    single_simb_undisp(i,:) = conv(single_simb_undisp(i,:), base_impulse_normalized, "same");
end

figure
plot(angle(ff))
grid on

figure
plot(real(tt))
grid on

figure
plot(abs(signal(3,:)))

figure(abs())

figure

