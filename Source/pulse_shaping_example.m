close all
clear variables

Nsym = 6;
beta = 0.5;
sampsPerSym = 20;

rctFilt = comm.RaisedCosineTransmitFilter(...
    Shape="Normal", ...
    RolloffFactor = beta, ...
    FilterSpanInSymbols=Nsym,...
    OutputSamplesPerSymbol=sampsPerSym);

b = coeffs(rctFilt);
rctFilt.Gain = 1/max(b.Numerator);

rctFilt.impz()
B = 40e9;
Fs = B * sampsPerSym; %Sampling Frequency
DataL = 2*Nsym + 1;
R = B;

hStr = RandStream('mt19937ar',Seed=0);
I = 2*randi(hStr,[0 1],DataL,1)-1;
Q = 2*randi(hStr,[0 1],DataL,1)-1;
%x = transpose([1, i, -i, 1, -1, -i,  i, 1, -i, i, -i, 1, -1]);
x = (I + 1i*Q);
tx = (1e9) * (0: DataL - 1) / R; % time vector in nanosecond

yo = rctFilt([x; zeros(Nsym/2,1)]);

%to = 1000000 * (0: (DataL+Nsym/2)*sampsPerSym - 1) / Fs;
to = (1e9) * (0: DataL*sampsPerSym - 1) / Fs;

fltDelay = Nsym / (2*R);
yo = yo(fltDelay*Fs+1:end);

fig1 = figure;
stem(tx, abs(x), 'kx'); hold on;

plot(to, abs(yo), 'b-');
plot(to, angle(yo), "b-")