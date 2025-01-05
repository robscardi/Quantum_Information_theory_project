close all
clear variables

%To modify the symbol in consideration: line 18
%To modify the local oscillator (laser) parameters: line 45
%To modify dispesion: line 130
%To modify constellation parameter (n bit): line 60
%To modify observation time: line 41

%% UNIT OF MEASURMENT
um = 1e-6;
nm = 1e-9;
pm = 1e-12;
MHz = 1e6;
kHz = 1e3;
km = 1e3;
ps = 1e-12;
mW = 1e-3;

%% CONSTANTS
symbol = (-1+17i);
B = 400e6;
c = 299792458;
h = 6.62607015e-34;

fiber_width = 9*um;

Aeff = pi*(fiber_width/2).^2;
eps0 = 8.8541878188e-12;
eps = (1.46).^2*eps0;

T = 1/B;
sps = 100;
ts = T/sps;
span = 16;

total_duration = span*T;
fs = 1/ts;
L = span*sps;
f = fs/L*(-L/2:L/2-1);
obs_time = 10*ts;
b = 0.2;
eta = 1;
maximum_field = 1e1; %Incoming signal field

%% LO DATA
lo.linewidth = 1*kHz;
lo.PSD = -40;
lo.lambda = 1550*nm;
lo.field = 1e3;
fc = c/lo.lambda;

phot_energy = h*fc;

sample_num = ceil(obs_time/ts);
sample_num = (mod(sample_num, 2) == 0)*(sample_num) + (mod(sample_num,2)==1)*(sample_num+1);


%% CONSTELLATION DATA
n_bit = 10;
M = 2^n_bit;
total_symbols = qammod(0:M-1, M);
symbol_vec = unique(real(total_symbols));
max_symbol = max(symbol_vec);

%% TEST AVG NUMBER OF PHOTON FOR 1+0i SYMBOL
% FOR UNCALIBRATED DETECTION
middle = sps*span/2 +1;
half_sample_num = ceil(sample_num/2);

base_impulse = 1*rcosdesign(b, span, sps);
base_impulse_normalized = base_impulse/(max(abs(fft(base_impulse))));

arrival_impulse = conv(base_impulse, base_impulse_normalized, "same");

%arrival_impulse_energy = trapz(abs(arrival_impulse).^2); %test symbol
%energy
symbol_sample = arrival_impulse(middle-half_sample_num: middle+half_sample_num);

%figure
%plot(arrival_impulse) 
avg_p = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field + lo.field).^2)*ts/phot_energy/obs_time;
avg_m = eta*eps*Aeff*trapz(abs(symbol_sample*maximum_field - lo.field).^2)*ts/phot_energy/obs_time;

max_phot_p = (avg_p-avg_m);

%% DISPERSION CALCULATION

lambda_vector = c./(f+fc);
Communication_lenght = 9*km;

% Calculation for optical fibre effective refractive index
% for this snippet to work, download
% https://it.mathworks.com/matlabcentral/fileexchange/27819-optical-fibre-toolbox
% and add to path
% materials = {'sm800core'; 'silica'};
%     fibre = struct(...
%     'materials', {materials});
% 
%     argument = struct(...
%     'type', 'wvl',... % calculate vs. wavelength
%     'harmonic', 1,... % required
%     'min', 1400,... % calculate from
%     'max', 2000 ...
%     ); % calculate to
% 
%     modeTask = struct(...
%     'nu', [0],... % first modal index
%     'type', {'te'},... % mode types
%     'maxmode', 1,... % how many modes of each type and NU to calculate
%     'diameter', 9);%,... % parameter, structure diameter, if argument is wavelength
%     %'region', 'cladding');
% 
%     infomode = false;
% 
%     modes = buildModes(argument, fibre, modeTask, infomode);
% 
%     neff = modes.NEFF;
%     l = modes.ARG;
% 
%     neff_interpolated = interp1(l, neff, lambda_vector/nm);
%     beta = zeros(size(f));


% for i = 1:length(lambda_vector)
%     neff_ = neff_interpolated(i);
%     beta(end+1-i) = 2*pi*neff_/lambda_vector(i);
% end
    
% MODIFY TO REGULATE DISPERSION
D = 2*(ps/(nm*km));
beta = D*((lo.lambda.*f).^2*pi/c);

ff = exp(-1i*beta*Communication_lenght);

tt = fftshift(fft(ff/length(ff)));

PPD = parallel.pool.PollableDataQueue;

%% SIMULATION

num_test = 100000;
tic;
parfor k =1:num_test
    random_indices = randi(length(total_symbols), 1, span+1);
    symbols = total_symbols(random_indices);
    symbols(span/2 +1 ) = symbol;
    [x, y]= optical_channel_func(lo,b, span, symbols, maximum_field, B, sps, sample_num, eta, tt, obs_time);
    send(PPD, x+1i*y);
end
toc

tic;
figure
hold on

inphase = zeros(1, num_test);
inquadrature = zeros(1, num_test);
syms = zeros(1, num_test);

for k = 1:num_test    
    z = poll(PPD);
    syms(k) = z;
    inphase(k) = real(z);
    inquadrature(k) = imag(z);
end

scatter(inphase, inquadrature, Color="blue");
toc

%% COMMENT: UNCALIBRATE DETECTION 
% TO SEE THE EFFECT OF DISPERSION WITHOUT THE CALIBRATION RUN, 
% USING ONLY THE EXPECTED MEAN (DOESN'T RISULT IN CORRECT VALUES FOR HIGH
% DISPERSION
 
inphase_freq = zeros(1,length(symbol_vec));
inquadr_freq = zeros(1,length(symbol_vec));


mean_inphase = mean(inphase);
mean_inquadrature = mean(inquadrature);

s = decode_qam(syms/max_phot_p, max_symbol);
    inphas = s(1,:);
    inquad = s(2,:);

for k = 1:length(symbol_vec)
    [~, index] = find(inphas == symbol_vec(k));
    inphase_freq(k) = length(index);

    [~, index] = find(inquad == symbol_vec(k));
    inquadr_freq(k) = length(index);
end

figure
title("Uncalibrated Distribution")
hold on
plot(symbol_vec, inphase_freq)
plot(symbol_vec, inquadr_freq)

figure 
title("Hard Decoding Uncalibrated")
hold on
scatter(inphas, inquad)

figure
hold on 
title("Uncalibrated Constellation")
scatter(real(total_symbols), imag(total_symbols), 'yellow', '*')
scatter(inphase/max_phot_p, inquadrature/max_phot_p, 'blue')


%% CALIBRATION RUN

symbol = 1+1i;

parfor k =1:num_test
    random_indices = randi(length(symbol_vec), 1, span+1);
    symbols = symbol_vec(random_indices);
    symbols(span/2 +1 ) = symbol;
    [x, y]= optical_channel_func(lo,b, span, symbols, maximum_field, B, sps, sample_num, eta, tt, obs_time);
    send(PPD, x+1i*y);
end

tinphase = zeros(1, num_test);
tinquadrature = zeros(1, num_test);
tsyms = zeros(1, num_test);

for k = 1:num_test    
    x = poll(PPD);
    tsyms(k) = x;
    tinphase(k) = real(x);
    tinquadrature(k) = imag(x);
end

t_mean_inphase = mean(tinphase);
t_mean_inquadrature = mean(tinquadrature);

s2 = decode_qam(inphase/t_mean_inphase+ ...
    1i*inquadrature/t_mean_inquadrature,max_symbol);
    
    t_inphase = s2(1,:);
    t_inquad = s2(2,:);
    
    t_inphase_freq = zeros(1,length(symbol_vec));
    t_inquadr_freq = zeros(1,length(symbol_vec));

for k = 1:length(symbol_vec)
    [~, index] = find(t_inphase == symbol_vec(k));
    t_inphase_freq(k) = length(index);

    [~, index] = find(t_inquad == symbol_vec(k));
    t_inquadr_freq(k) = length(index);
end

% figure
% hold on 
% title("Calibrated constellation: ")
% scatter(real(total_symbols), imag(total_symbols), 'yellow', '*')
% scatter(tinphase/t_mean_inphase, tinquadrature/t_mean_inquadrature, 'blue')

figure
hold on
title("Calibrated distribution")
plot(symbol_vec, t_inphase_freq, DisplayName="Inphase")
plot(symbol_vec, t_inquadr_freq, DisplayName="Inquadr")
legend on

figure
hold on
title("Calibrated constellation : " + 2^n_bit + "QAM")
plot(real(total_symbols), imag(total_symbols), '*', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
scatter(inphase/t_mean_inphase, inquadrature/t_mean_inquadrature, 'blue')
grid on


function [s] = decode_qam(symbol, max_symbol)
    p = real(symbol);
    q = imag(symbol);
    m = max_symbol;
    s(1,:) = nearest_odd(p);
    s(2,:) = nearest_odd(q);
    s(s>m) = m;
    s(s<-m) = -m;
end