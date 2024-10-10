%% ORIGINAL OPTICAL CHANNEL WITH NORMAL RAISED_COSINE
close all
clear variables

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
B = 25e9;           %Symbol rate
c = 299792458;

D = 17*ps/(km*nm);
h = 6.62607015e-34;


maximum_field=1e5; %V/m

%% LO DATA
lo.linewidth = 1*MHz;
lo.PSD = -80;
lo.lambda = 1550*nm; 
fc = c/lo.lambda; 

b2 = -((lo.lambda^2) /(2*pi*c))*D;

fiber_width = 9*um;

Aeff = pi*(fiber_width/2).^2;
eps0 = 8.8541878188e-12;
eps = (1.46).^2*eps0;


%% PULSE SHAPING
% ts = 0.01;
% t = (-3:ts:3); % Time vector
% beta = 0.5;   % Roll-off factor
% T = 1;         % Symbol period (for example, T = 1)
% 
% sig = zeros(7, length(t));
% hh = zeros(1, length(t));
% 
% for i = 0:1:6
%     sig(i+1,:) = raised_cosine_impulse(t+i-3, beta, T)';
%     hh = hh + sig(i+1,:);
% end

T = 1/B;
beta = 0.5;
ts = T/1000;
span = 8;
sps = 100;
base_impulse = rcosdesign(beta, span, sps);

symbols = ones(span);

t = linspace(-((span+1)/2)*T,  ((span+1)/2)*T, sps*span +1);
base_impulse = base_impulse/max(base_impulse);
signal = zeros(span+1, length(t));
total_sign = zeros(1, length(t));

for i = 0:(span)
    total_sign(i*sps + 1) = 1;
    signal(i+1, i*sps +1 ) = 1;
    signal(i+1, :) = conv(base_impulse, signal(i+1, :), "same");
end

total_sign = conv(base_impulse, total_sign, "same");

figure
hold on
plot(t, total_sign)
plot(t, signal)



%% DISPERSION SIMULATION

Communication_length = 10*km;


dT = b2*Communication_length*2*B;

disp_h = zeros(7, length(t));
disp_hh = zeros(1, length(t));

a = rand(1, 7);

for i = 0:1:6
    disp_h(i+1,:) = a(i+1)*raised_cosine_impulse(t+i-3, beta, T+dT*B)';
    disp_hh = disp_hh + disp_h(i+1,:);
end

fig8 = figure();
hold on
plot(t/B, disp_h);
plot(t/B, disp_hh, LineWidth=1.2);
xlabel('Time');
ylabel('Amplitude');
title('Raised Cosine Impulse Response');
grid on;


% 
%     Communication_lenght = 1000*km;
% 
%     materials = {'sm800core'; 'silica'};
%     fibre = struct(...
%     'materials', {materials});
% 
%     argument = struct(...
%     'type', 'wvl',... % calculate vs. wavelength
%     'harmonic', 1,... % required
%     'min', minlambda*1e6,... % calculate from
%     'max', maxlambda*1e6); % calculate to
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
%     fig4 = figure();
%     showModes(modes, 'Modal Dispersion in Fiber')
% 
%     neff = modes.NEFF;
%     l = modes.ARG;
% 
%     neff_interpolated = interp1(l, neff, lambda_vector/um);
%     beta = zeros(size(f));
% 
% 
%     for i = 1:length(lambda_vector)
%         neff_ = neff_interpolated(i);
%         beta(end+1-i) = 2*pi*neff_/lambda_vector(i);
%     end
% 
%     %hh_f = hh_f .* exp(1i.*beta.*Communication_lenght);
%     dt = exp(1i*beta*Communication_lenght);
% 
% 
%     fig7 = figure();
%     %plot(f, abs(hh_f))
%     dispersed_signal = ifft(hh_f);  
%     fig8 = figure();
%     hold on
%     plot(t,dispersed_signal);
%     plot(t,hh);
% 
% 
% 
% end


%% PHASE NOISE
cycle_time = 1/fc;
sampling_time = 50*1000; %in samples


fmin = fc - lo.linewidth/2;
fmax = fc + lo.linewidth/2;

D = 10^(lo.PSD/10);
var = D*((1/fmin)-(1/fmax));
std = sqrt(var);
fs = 3*fc;
ts = 10000/fs;

t = ts*(0:(sampling_time-1));

phi = zeros(sampling_time, 1);
random_walk = std*randn(1, length(phi)-1);

for i = 1:1:length(phi)-1
    phi(i+1) = phi(i) + random_walk(i);
end

lo.wave = sin(2*pi*fc*t + phi');
y = fftshift(fft(lo.wave))/sampling_time;
fig3 = figure();
semilogy(fs/sampling_time*(-sampling_time/2 : (sampling_time/2 -1)), abs(y));


%% SYMBOL CALC


t = t- t(ceil(length(t)/2));
input_phase= 2*pi*rand();
input_power= maximum_field*rand();

lo.field = maximum_field;

phot_energy = h*fc;

maximum_power = (maximum_field*fiber_width.^2) * pi/2 %Gaussian mode

max_avg_numb_photon = maximum_power/phot_energy

a = rand(1, 7);
phi_simbol = 2*pi*rand(1,7);

disp_h = zeros(7, length(t));
disp_hh = zeros(1, length(t));

for i = 0:1:6
    disp_h(i+1,:) = a(i+1).*exp(1i*phi_simbol(i+1)).*raised_cosine_impulse(t+i-3, beta, T+dT*B)';
    disp_hh = disp_hh + disp_h(i+1,:);
end
 
fig10 = figure();
hold on
plot(abs(disp_hh), t);
plot(angle(disp_hh), t);

symbol_amp = rand();
symbol_phase = 2*pi*rand();

%symbol_in_phase = po*maximum_field*cos(input_phase)/phot_energy;
%symbol_in_quadrature = po*maximum_field*sin(input_phase)/phot_energy;


phot_avg_in_phase =      ((maximum_field*cos(symbol_phase)*symbol_amp).^2)*(2*eps*Aeff)/(phot_energy);
phot_avg_in_quadrature = ((maximum_field*sin(symbol_phase)*symbol_amp).^2)*(2*eps*Aeff)/(phot_energy);



constellation = figure();

plot(phot_avg_in_phase, phot_avg_in_quadrature, Marker="*", Color="yellow")

%plot(symbol_in_phase/max_avg_numb_photon, symbol_in_quadrature/max_avg_numb_photon,"Color","yellow", "Marker",".", "MarkerSize", 10);
hold on

%plot(I_avg_num_of_phot/max_avg_numb_photon, Q_avg_num_of_phot/max_avg_numb_photon,"Color","green", "Marker",".", "MarkerSize", 10);




for j = 1:1000

    a = ones(1,7)*symbol_amp; %rand(1, 7);
    phi_simbol = symbol_phase*ones(1,7); % rand(1, 7);
    a(3) = symbol_amp;
    phi_symbol(3) = symbol_phase;

    disp_h = zeros(7, length(t));
    disp_hh = zeros(1, length(t));

    for i = 0:1:6
        disp_h(i+1,:) = exp(1i*phi_simbol(i+1))*a(i+1)*raised_cosine_impulse(t+i-3, beta, T+dT*B)';
        disp_hh = disp_hh + disp_h(i+1,:);
    end
    %% IN PHASE

    lo.laser = lo.field*exp(1i*phi);
    avg_p = eps*Aeff*trapz(abs(disp_hh*maximum_field + lo.laser').^2)*ts/phot_energy;
    avg_m = eps*Aeff*trapz(abs(disp_hh*maximum_field - lo.laser').^2)*ts/phot_energy;

    k_p = poissrnd(avg_p);
    k_m = poissrnd(avg_m);

    I_phot = k_p - k_m;
   
    %% IN QUADRATURE
    
    lo.laser = lo.field*exp(1i*(phi + pi/2));
    avg_p = eps*Aeff*trapz(abs(disp_hh*maximum_field + lo.laser').^2)*ts/phot_energy;
    avg_m = eps*Aeff*trapz(abs(disp_hh*maximum_field - lo.laser').^2)*ts/phot_energy;

    k_p = poissrnd(avg_p);
    k_m = poissrnd(avg_m);

    Q_phot = k_p - k_m;
    
    plot(I_phot, Q_phot,"Color","blue", "Marker",".", "MarkerSize", 10);

end

xlabel("In Phase")
ylabel("In Quadrature")

th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);
%plot(xunit, yunit);

grid on
