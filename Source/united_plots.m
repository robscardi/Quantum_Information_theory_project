clear variables
close all

um = 1e-6;
nm = 1e-9;
pm = 1e-12;
MHz = 1e6;
kHz = 1e3;
km = 1e3;
ps = 1e-12;
mW = 1e-3;
x = linspace(1, 5e2, 30);
b = 0.2;
eta = 1;
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
lambda = 1550*nm;
fc = c/lambda;

phot_energy = h*fc;
sample_num = ceil(obs_time/ts);
sample_num = (mod(sample_num, 2) == 0)*(sample_num) + (mod(sample_num,2)==1)*(sample_num+1);



middle = sps*span/2 +1;
half_sample_num = ceil(sample_num/2);

base_impulse = 1*rcosdesign(b, span, sps);
base_impulse_normalized = base_impulse/(max(abs(fft(base_impulse))));

arrival_impulse = conv(base_impulse, base_impulse_normalized, "same");

%arrival_impulse_energy = trapz(abs(arrival_impulse).^2);
symbol_sample = arrival_impulse(middle-half_sample_num: middle+half_sample_num);


%figure
%plot(arrival_impulse)
k = 1;
avg_p = zeros('like', x);
avg_m = zeros('like', x);
for g= x
    avg_p(k) = eta*eps*Aeff*trapz(abs(symbol_sample*10 + g).^2)*ts/phot_energy/obs_time;
    avg_m(k) = eta*eps*Aeff*trapz(abs(symbol_sample*10 - g).^2)*ts/phot_energy/obs_time;
    k = 1+k;
end
max_phot = (avg_p-avg_m);

g = zeros(10, length(x));
q = zeros(10, length(x));
data_string = "17_disp_100_km_1khz_-80dbm";

for i=2:8
    a = load("..\Data\" + data_string + "\"+ i +"_bit.mat");
    g(i,:) = a(1).MI_vector;
    a = load("..\Data\" + data_string + "\"+ "QBER_"+ i +"_bit.mat");
    q(i,:) = a(1).QBER;
end

C_h = ((max_phot+1).*log2(max_phot+1)- (max_phot).*log2(max_phot))*0.5;
C_s2 = log2(1+2*max_phot)*0.5;

figure("Name", "Mutual Information")
hold on
for i=2:10
    plot(max_phot, g(i,:), "DisplayName",2^i +"QAM", LineWidth=2, Color="#0072BD")
end

for i=2:9
     text(13 + (-1)^i, g(i,end)-0.05, 2^i +"QAM", fontsize=15)
end
i = 10;
    text(13 + (-1)^i, g(i,end)+0.04, 2^i +"QAM", fontsize=15)

plot(max_phot, C_h, DisplayName="Holevo Bound", LineWidth=2)
text(12, C_h(end)-0.05, "Holevo Bound", fontsize=25)
plot(max_phot, C_s2, DisplayName="Shannon Bound", LineWidth=2)
text(12, C_s2(end)-0.38,"Shannon Bound", fontsize=25)
%legend("FontSize",15)
xlim([2, max_phot(end)])
xlabel("# Photon", FontSize=45)
ylabel("Mutual Information [bits]", FontSize=45)
grid on
ax = gca; % Get the current axis
ax.FontSize = 30;

figure("Name", "Bit error rate")
for i=2:10
    semilogy(max_phot, q(i,:), "DisplayName",2^i +"QAM", LineWidth=2)
    hold on
end
xlim([2, max_phot(end)])
xlabel("# Photon", FontSize=45)
ylabel("Bit Error Rate", FontSize=45)
legend("FontSize",25, 'Location', 'best', 'NumColumns', ceil(numel(findall(gca, 'Type', 'line')) / 2));
grid on
ax = gca; % Get the current axis
ax.FontSize = 30;