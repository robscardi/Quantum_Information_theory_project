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
maximum_field = 1e1;

%% LO DATA
lo.linewidth = 1*kHz;
lo.PSD = 20;
lo.lambda = 1550*nm;
lo.field = 1e1;
fc = c/lo.lambda;

phot_energy = h*fc;

sample_num = ceil(obs_time/ts);
sample_num = (mod(sample_num, 2) == 0)*(sample_num) + (mod(sample_num,2)==1)*(sample_num+1);

N_e_vector = 30;
E_vector = linspace(1e1, 5e2, N_e_vector);
MI_vector = zeros(1,N_e_vector);

%% CONSTELLATION DATA
n_bit = 10;
M = 2^n_bit;
total_symbols = qammod(0:M-1, M);
symbol_vec = unique(real(total_symbols));
max_symbol = max(symbol_vec);
num_test = 1000;

%% DISPERSION CALCULATION

lambda_vector = c./(f+fc);
Communication_lenght = 9*km;
    
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

    D = 0*(ps/(nm*km));
    beta = D*((lo.lambda.*f).^2*pi/c);

    ff = exp(-1i*beta*Communication_lenght);
    
    tt = ifftshift(ifft(ifftshift(ff)));

PPD = parallel.pool.PollableDataQueue;
QBER = zeros(1, length(E_vector));
for g =1:N_e_vector
    lo.field = E_vector(g);
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

    cross_prob_matrix_p = zeros(length(symbol_vec));
    cross_prob_matrix_q = zeros(length(symbol_vec));

    %% SIMULATION


    for symbol = total_symbols

        parfor k =1:num_test
            random_indices = randi(length(total_symbols), 1, span+1);
            symbols = total_symbols(random_indices);
            symbols(span/2 +1 ) = symbol;
            [x, y]= optical_channel_func(lo,b, span, symbols, maximum_field, B, sps, sample_num, eta, tt, obs_time);
            send(PPD, x+1i*y);
        end

        inphase = zeros(1, num_test);
        inquadrature = zeros(1, num_test);
        syms = zeros(1, num_test);

        for k = 1:num_test
            z = poll(PPD);
            syms(k) = z;
            inphase(k) = real(z);
            inquadrature(k) = imag(z);
        end

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

        [~,col_p] = find(symbol_vec == real(symbol));
        [~,col_q] = find(symbol_vec == imag(symbol));
        
        QBER_P = (num_test - t_inphase_freq(col_p));
        QBER_Q = (num_test - t_inquadr_freq(col_q));
        
        QBER(g) = QBER(g) + QBER_P*0.5 + QBER_Q*0.5;


        cross_prob_matrix_p(col_p, :) = cross_prob_matrix_p(col_p, :) + t_inphase_freq;
        cross_prob_matrix_q(col_q, :) = cross_prob_matrix_q(col_q, :) + t_inquadr_freq;
    end

    %% NORMALIZE PROBABILITIES

    for i = 1:length(symbol_vec)
        N = sum(cross_prob_matrix_p(i, :));
        if N > 0
            cross_prob_matrix_p(i, :) = cross_prob_matrix_p(i, :)/N;
        end
        N = sum(cross_prob_matrix_q(i, :));
        if N > 0
            cross_prob_matrix_q(i, :) = cross_prob_matrix_q(i, :)/N;
        end
    end

    %% CALCULATE MUTUAL INFORMATION


    v = 1/(2*max_symbol);%/2;
    X_p = zeros(1, length(symbol_vec));
    
    for k = 1:length(symbol_vec)
        for l = symbol_vec
            X_p(k) = X_p(k) + exp(-v*(l^2 + symbol_vec(k)^2));
        end
    end
    N = sum(X_p);

    X_p = X_p/N; 

    MI_Phase = mutual_information(X_p, cross_prob_matrix_p);
    MI_Quadr = mutual_information(X_p, cross_prob_matrix_q);

    MI_vector(g) = (MI_Phase*0.5 + MI_Quadr*0.5)*0.5;
    
end 

QBER = QBER/(num_test*length(total_symbols));

figure
hold on
plot(E_vector, MI_vector)

figure
hold on
plot(E_vector, QBER);

function [s] = decode_qam(symbol, max_symbol)
    p = real(symbol);
    q = imag(symbol);
    m = max_symbol;
    s(1,:) = nearest_odd(p);
    s(2,:) = nearest_odd(q);
    s(s>m) = m;
    s(s<-m) = -m;
end