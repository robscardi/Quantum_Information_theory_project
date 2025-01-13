function [result] = fieldPower(F)
% Power in the mode field.

%% Create mesh grid 
FG = fieldGrid(F);

% In region 1
CrossProduct1 = cross (F.E1, conj(F.H1), 3) + cross (conj(F.E1), F.H1, 3); 
% Power calculated with this equation, taken from Gubsky 2005, in test on
% 2009-12-01-KK gave the same result as with equation 0.5ExH* as in Yariv.

% Calculating the discretization steps for R and PHI
%  dr = max(max(R)) ./ size (R, 2);
%  dphi = 2*pi ./ size (PHI, 1);

% The surface elements in dStot are calculated for all r 
% dS = (dr * dphi) .* R; 

% ERROR: ZProduct = CrossProduct(:,:,3) .* dS;
value1 = 0.25 * sum(sum(CrossProduct1(:,:,3) .* FG.ds1));

% In region 2
value2 = 0;
if ~isempty(F.E2)
    CrossProduct2 = cross (F.E2, conj(F.H2), 3) + cross (conj(F.E2), F.H2, 3);
    value2 = 0.25 * sum(sum(CrossProduct2(:,:,3) .* FG.ds2));
end;

result = value1 + value2;

% Generating a cylinder - it is used only for plotting the fiber 
% A = zeros(1,20);
%       for i = 1:20
%           A(i) = radius;
%       end;
% [Cil_x, Cil_y, Cil_z] = cylinder(A);
% 
%  figure; 
%  % Drawing the cylinder
%  surf(Cil_x, Cil_y, max(max(ZProduct)) * Cil_z, 'EdgeAlpha',1,'FaceAlpha',0);
%  hold on; 
%  % Drawing the product (F_j x G_j^* + F_j^* x G_j)z in each point (eq. 2, Grubsky & Savchenko)
%  surf( X, Y, ZProduct, 'EdgeAlpha',0,'FaceAlpha',1);
%  title('Scalar Product');xlabel ('X'); ylabel('Y'); axis on; view([0 90]);
%  colorbar;
% 
