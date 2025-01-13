function [code] = displayField2(F, d)
% Creates plots of the mode field for two-layer structure

assert(nargin == 2, 'Wrong number of input arguments\n', mfilename);

% Correction for nice pictures (no white sector)
F.FG.PHI1 = cat(1, F.FG.PHI1, F.FG.PHI1(1, :) + 2*pi); 
F.FG.R1 = cat(1, F.FG.R1, F.FG.R1(1, :)); 
F.E1 = cat(1, F.E1, F.E1(1, :, :)); 
F.H1 = cat(1, F.H1, F.H1(1, :, :)); 

F.FG.PHI2 = cat(1, F.FG.PHI2, F.FG.PHI2(1, :) + 2*pi); 
F.FG.R2 = cat(1, F.FG.R2, F.FG.R2(1, :)); 
F.E2 = cat(1, F.E2, F.E2(1, :, :)); 
F.H2 = cat(1, F.H2, F.H2(1, :, :)); 

% Convert to Cartesian coordinates
[X1,Y1] = pol2cart(F.FG.PHI1,F.FG.R1);
[X2,Y2] = pol2cart(F.FG.PHI2,F.FG.R2);

EX1 = F.E1(:,:,1) .* cos(F.FG.PHI1) - F.E1(:,:,2) .* sin(F.FG.PHI1);
EY1 = F.E1(:,:,1) .* sin(F.FG.PHI1) + F.E1(:,:,2) .* cos(F.FG.PHI1);
EZ1 = F.E1(:,:,3);
EX2 = F.E2(:,:,1) .* cos(F.FG.PHI2) - F.E2(:,:,2) .* sin(F.FG.PHI2); 
EY2 = F.E2(:,:,1) .* sin(F.FG.PHI2) + F.E2(:,:,2) .* cos(F.FG.PHI2); 
EZ2 = F.E2(:,:,3);

% Calculating the total value of the field at each coordinate
E1 = sqrt(abs(EX1).^2 + abs(EY1).^2 + abs(EZ1).^2);
E2 = sqrt(abs(EX2).^2 + abs(EY2).^2 + abs(EZ2).^2);

% % Calculating the intensity at each coordinate. Note: for
% % electric fields, this gives a value that is proportional, but not equal,
% % to the optical intensity. 
% I = E.^2;

% Make sure that the values to be plotted are real. If not, each
% component is made equal to its absolute value. 
if ~isreal(EY1)
    EY1 = real(EY1 * 1i);
end;  

if ~isreal(EY2)
    EY2 = real(EY2 * 1i);
end;  

if ~isreal(EX1)
    EX1 = real(EX1 * 1i);
end;  

if ~isreal(EX2)
    EX2 = real(EX2 * 1i);
end;  

if ~isreal(EZ1)
    EZ1 = real(EZ1 * 1i);
end;  

if ~isreal(EZ2)
    EZ2 = real(EZ2 * 1i);
end;  


%%

figure;

subplot(2,2,1);
hold on
surf( X1, Y1, EX1,'EdgeAlpha',0,'FaceAlpha',1);
surf( X2, Y2, EX2,'EdgeAlpha',0,'FaceAlpha',1);
theta = 0:.1:2*pi;
polar(theta,ones(size(theta)) * d/2)
title('E_x (a.u.)');xlabel ('X (\mum)'); ylabel('Y(\mum)'); axis on; view([0 90]);
colorbar;

subplot(2,2,2);
hold on
surf( X1, Y1, EY1, 'EdgeAlpha',0,'FaceAlpha',1);
surf( X2, Y2, EY2,'EdgeAlpha',0,'FaceAlpha',1);
polar(theta,ones(size(theta)) * d/2)
title('E_y (a.u.)');xlabel ('X (\mum)'); ylabel('Y(\mum)'); axis on; view([0 90]);
colorbar;

subplot(2,2,3);
hold on
surf( X1, Y1, EZ1, 'EdgeAlpha',0,'FaceAlpha',1);
surf( X2, Y2, EZ2,'EdgeAlpha',0,'FaceAlpha',1);
polar(theta,ones(size(theta)) * d/2)
title('E_z (a.u.)');xlabel ('X (\mum)'); ylabel('Y(\mum)'); axis on; view([0 90]);
colorbar;

subplot(2,2,4);
hold on
surf( X1, Y1, E1, 'EdgeAlpha',0,'FaceAlpha',1);
surf( X2, Y2, E2, 'EdgeAlpha',0,'FaceAlpha',1);
polar(theta,ones(size(theta)) * d/2)
title('E_{total} --- absolute value');xlabel ('X (\mum)'); ylabel('Y(\mum)'); axis on; view([0 90]);
colorbar;

%% Display a 3D image of "intenisty distribution"
% Calculate the intensity at each coordinate. Note: for
% electric fields, this gives a value that is proportional, but not equal,
% to the optical intensity. 
E1 = sqrt(abs(EX1).^2 + abs(EY1).^2 + abs(EZ1).^2);
I1 = E1.^2;
E2 = sqrt(abs(EX2).^2 + abs(EY2).^2 + abs(EZ2).^2);
I2 = E2.^2;

% Create a cylinder having the same radius as the fiber
A = zeros(1,20);
      for i = 1:20
          A(i) = d/2;
      end;
[Cil_x, Cil_y, Cil_z] = cylinder(A);
 
figure; 
% Plot the cylinder
height = max([max(max(I1)) max(max(I2))]);
surf(Cil_x, Cil_y, height * Cil_z, 'EdgeAlpha',1,'FaceAlpha',0);
hold on;  
surf(Cil_x, Cil_y, height * Cil_z, 'EdgeAlpha',1,'FaceAlpha',0);
% plot the "intensity" 
surf(X1,Y1, I1, 'EdgeAlpha',0,'FaceAlpha',1);
surf(X2,Y2, I2, 'EdgeAlpha',0,'FaceAlpha',1);
title('Intensity (a.u.)');xlabel ('X (\mum)'); ylabel('Y(\mum)'); axis on; view([-30 60]);
%colorbar;

% test2bis=figure;
% quiver(X,Y, EX, EY,0.5);

% test3bis=figure;
% Z_tot = zeros(size(X));
% quiver3(X,Y, Z_tot, EX, EY, abs(EZ), 0.5);

% quiver(X, Y, EX, EY, 0.5);
% hold on;
% phi = 0:0.1:pi;
% plot(a * cos(phi), a * sin(phi));
% plot(a * cos(phi), -a * sin(phi));
% return;
