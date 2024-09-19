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


symbol = (1+1i)/(sqrt(2));
B = 25e9;
c = 299792458;

T = 1/B;
sps = 200;
ts = T/sps;
span = 16;

total_duration = span*T;
fs = 1/ts;
L = span*sps;
f = fs/L*(-L/2:L/2-1);
obs_time = 5*ps; 

%% LO DATA
lo.linewidth = 1*MHz;
lo.PSD = -80;
lo.lambda = 1550*nm;
lo.field = 1e9;
fc = c/lo.lambda;

lambda_vector = c./(f+fc);

Communication_lenght = 1000*km;

    materials = {'sm800core'; 'silica'};
    fibre = struct(...
    'materials', {materials});

    argument = struct(...
    'type', 'wvl',... % calculate vs. wavelength
    'harmonic', 1,... % required
    'min', 1400,... % calculate from
    'max', 2000 ...
    ); % calculate to

    modeTask = struct(...
    'nu', [0],... % first modal index
    'type', {'te'},... % mode types
    'maxmode', 1,... % how many modes of each type and NU to calculate
    'diameter', 9);%,... % parameter, structure diameter, if argument is wavelength
    %'region', 'cladding');

    infomode = false;

    modes = buildModes(argument, fibre, modeTask, infomode);

    neff = modes.NEFF;
    l = modes.ARG;

    neff_interpolated = interp1(l, neff, lambda_vector/nm);
    beta = zeros(size(f));


    for i = 1:length(lambda_vector)
        neff_ = neff_interpolated(i);
        beta(end+1-i) = 2*pi*neff_/lambda_vector(i);
    end

    ff = exp(1i*beta*Communication_lenght);
    
    tt = ifft(ff)/L;

D = parallel.pool.PollableDataQueue;

num_test = 1000;
tic;
for k =1:num_test
    [x, y]= optical_channel_func(lo, 16, symbol, 1e8, 25e9, sps, 1, tt, obs_time);
    send(D, x+1i*y);
end
toc


tic;
figure
hold on

inphase = nan(1, num_test);
inquadrature = nan(1, num_test);
for k = 1:num_test    
    x = poll(D);
    inphase(k) = real(x);
    inquadrature(k) = imag(x);
end

scatter(inphase, inquadrature, Color="blue");

min_inphase = min(inphase);
max_inphase = max(inphase);
inphase_freq = nan(1, max_inphase - min_inphase+1);

min_inquad = min(inquadrature);
max_inquad = max(inquadrature);
inquadrature_freq = nan(1, max_inquad-min_inquad +1 );

for k = 1:num_test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if isnan(inphase_freq(inphase(k)-min_inphase+1)) 
        inphase_freq(inphase(k)-min_inphase+1) = 1;
    else 
        inphase_freq(inphase(k)-min_inphase+1) = inphase_freq(inphase(k)-min_inphase+1)+1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isnan(inquadrature_freq(inquadrature(k)-min_inquad+1)) 
        inquadrature_freq(inquadrature(k)-min_inquad+1) = 1;
    else 
        inquadrature_freq(inquadrature(k)-min_inquad+1) = inquadrature_freq(inquadrature(k)-min_inquad+1)+1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

figure
hold on
plot(min_inphase:max_inphase, inphase_freq)
plot(min_inquad:max_inquad, inquadrature_freq)

toc