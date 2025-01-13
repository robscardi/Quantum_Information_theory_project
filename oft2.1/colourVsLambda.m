function c = colourVsLambda(x)
% Returns a colour as corresponding to the given wavelength. 
% All parameters are set by Kotya and do not correspond to any
% agreements or sources.
%
% Wavelength is in nanometers. Returns ColorSpec [R G B] 
%
% (by) Karapetyan, 2008-2009.
% kotya.karapetyan@gmail.com
% http://agmeschede.iap.uni-bonn.de/

IRmax = 0.3;
Rmax = 1;
Gmax = 1;
Bmax = 1;
UVRmax = 0.5;
UVBmax = 0.4;

V = 200;
B = 400;
C = 450;
G = 550;
Y = 590;
R = 650;
IR = 800;

if x > IR
    c = [IRmax 0 0];
elseif x > R
    c = [Rmax - (x-R)*(Rmax-IRmax)/(IR-R) 0 0];
elseif x > Y
    c = [Rmax Gmax - (x-Y)*Gmax/(R-Y) 0];
elseif x > G
    c = [(x - G)*Rmax/(Y-G) Gmax 0];
elseif x > C
    c = [0 Gmax Bmax - (x-C)*Bmax/(G-C)];
elseif x > B
    c = [0 (x - B)*Gmax/(C-B) Bmax];
elseif x > V
    c = [UVRmax - (x-V)*UVRmax/(B-V) 0 UVBmax + (x - V)*(Bmax-UVBmax)/(B-V)];
else
    c = [UVRmax 0 UVBmax];
end;

for i=1:3
    if c(i)>1
        c(i)=1;
    end;
    if c(i)<0
        c(i)=0;
    end;
end;
assert(sum(c>=0) == 3 && sum(c<=1) == 3);

% UV = 0.300;
% UB = 0.400;
% BG = 0.500;
% GR = 0.600;
% IR = 0.800;
% 
% if lambda < UV
%     result = [0.05 0 0.05];
% elseif lambda < UB
%     result = [0.1 0 1];
% elseif lambda < BG
%     result = [0 0 1];
% elseif lambda < GR
%     result = [0 1 0];
% elseif lambda < IR
%     result = [1 0 0];
% else
%     result = [0.3 0 0];
% end;
% 
% 
% if lambda < 0.4
%     result = [0 0 1];
% elseif lambda < 0.7
%     result = [0 1 0];
% elseif lambda < 1.1
%     result = [1 0 0];
% else
%     result = [0 0 0];
% end;

