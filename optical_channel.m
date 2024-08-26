close all
clear variables

%% UNIT OF MEASURMENT
um = 1e-6;
nm = 1e-9;
MHz = 1e6;
kHz = 1e3;
km = 1e3;

B = 50e9;           %Symbol rate
c = 299792458;

%% LO DATA
lo.linewidth = 1*MHz;
lo.PSD = -80;
lo.lambda = 1550*nm; 
fc = c/lo.lambda; 

%% PULSE SHAPING
ts = 0.01;
t = (-3:ts:3); % Time vector
beta = 0.5;   % Roll-off factor
T = 1;         % Symbol period (for example, T = 1)

h = zeros(7, length(t));
hh = zeros(1, length(t));
for i = 0:1:6
    h(i+1,:) = raised_cosine_impulse(t+i-3, beta, T)';
    hh = hh + h(i+1,:);
end
fig1 = figure();
hold on
plot(t/B, h);
plot(t/B, hh);
xlabel('Time');
ylabel('Amplitude');
title('Raised Cosine Impulse Response');
grid on;

hh_f = fft(hh);

fs = B/(ts);
L = length(t);
f = fs/L*(-L/2:(L/2)-1);

fig5 = figure();
plot( f, abs(fftshift(hh_f))/L)


%% DISPERSION SIMULATION
dispersion = 1;
f = f + fc;

maxf = f(end);
minf = f(1);

maxlambda = c/minf;
minlambda = c/maxf;

lambda_vector = c./f;

fig6 = figure();
plot(lambda_vector, abs(hh_f))

if dispersion
    
    Communication_lenght = 1000*km;

    materials = {'sm800core'; 'silica'};
    fibre = struct(...
    'materials', {materials});
    
    argument = struct(...
    'type', 'wvl',... % calculate vs. wavelength
    'harmonic', 1,... % required
    'min', minlambda*1e6,... % calculate from
    'max', maxlambda*1e6); % calculate to
    
    modeTask = struct(...
    'nu', [0],... % first modal index
    'type', {'te'},... % mode types
    'maxmode', 1,... % how many modes of each type and NU to calculate
    'diameter', 9);%,... % parameter, structure diameter, if argument is wavelength
    %'region', 'cladding');

    infomode = false;
    
    modes = buildModes(argument, fibre, modeTask, infomode);
    fig4 = figure();
    showModes(modes, 'Modal Dispersion in Fiber')
    
    neff = modes.NEFF;
    l = modes.ARG;

    neff_interpolated = interp1(l, neff, lambda_vector/um);
    beta = zeros(size(f));


    for i = 1:length(lambda_vector)
        neff_ = neff_interpolated(i);
        beta(end+1-i) = 2*pi*neff_/lambda_vector(i);
    end

    %hh_f = hh_f .* exp(1i.*beta.*Communication_lenght);
    dt = exp(1i*beta*Communication_lenght);
    
    
    fig7 = figure();
    %plot(f, abs(hh_f))
    dispersed_signal = ifft(hh_f);  
    fig8 = figure();
    hold on
    plot(t,dispersed_signal);
    plot(t,hh);
    


end



%% PHASE NOISE
cycle_time = 1/fc;
sampling_time = 50*100; %in samples


fmin = fc - lo.linewidth/2;
fmax = fc + lo.linewidth/2;

D = 10^(lo.PSD/10);
var = D*((1/fmin)-(1/fmax));
std = sqrt(var);
fs = 3*fc;
ts = 1/fs;

t = ts*(0:(sampling_time-1));

phi = zeros(sampling_time, 1);
rand = std*randn(length(phi)-1);
phi(1) = 0;
for i = 1:1:length(phi)-1
    phi(i+1) = phi(i) + rand(i);
end

lo.wave = sin(2*pi*fc*t + phi');
fig2 = figure();
plot(t, lo.wave);
y = fftshift(fft(lo.wave))/sampling_time;
fig3 = figure();
semilogy(fs/sampling_time*(-sampling_time/2 : (sampling_time/2 -1)), abs(y));

I_comp = cos(phi);
Q_comp = sin(phi);



