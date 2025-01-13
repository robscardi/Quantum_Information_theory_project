%% Effective refractive index of three-layer modes 
% This tutorial demonstrates the use of the Optical Fibre Toolbox for
% calculation of the three-layer (core-cladding-surrounding) fibres. First,
% the specification of such structures is explained. Then the effective
% refractive index of the guided modes is calculated vs. fibre diameter
% (typical task for tapered microfibres).   
%
% (cc-by) K. Karapetyan et al., AG Meschede, Uni Bonn, 2011
%
% http://agmeschede.iap.uni-bonn.de | kotya.karapetyan@gmail.com

%% 
calculatemodes = true; % decide to enable or disable (for speed) calculation
if ~calculatemodes % if calculation is disabled the data should be previously stored in a mat file
    load tutorial3ls.mat modesMonerieCore modesTsaoTClad modesTwoClad modesTwoCore modesErdoganClad
end

%% Three layer fibre
% Usually, optical fibres are considered consisting of just two layers: the
% core and the cladding. Light field is well confined in the core and in
% the evanescent field around it in the cladding and does not reach the
% external surface of the cladding.
%
% However, there are cases when this model is not valid. One of them is the
% tapered optical microfibre, consisting of the untapered ends, tapers, and
% micrometer-diameter waist. In such fibres, the light is first guided as 
% usual, by the core in the untapered part. In the taper region, the light
% field expands and is no longer confined close to the core. At some point,
% the guidance is actually not by the core-cladding interface, but by the
% cladding-surrounding (usually --- air, liquids, or special coating having
% the refractive index lower than that of the cladding). In the waist the
% core is so small in diameter that it can be neglected in most cases and
% the two layer model can again be used, this time for
% cladding-surrounding. So, it is the taper region where the two-layer
% model is not enough to simulate light propagation in the fibre, and OFT
% provides the three-layer solutions for these cases.
%
%% Known solutions
% Simulation of the guided modes in the three-layer structure is an
% analytically challenging task. Several attempts have been made to publish
% the solutions in the literature. In 1976, BELANOV et al. have made the
% first known to me attempt to write the equations for the three-layer
% structure. However they have only solved them to the end for the weakly
% guiding fibre (linearly polarised modes approximation). It was also
% independently done by MONERIE in 1982. In 1989, TSAO has made an attempt
% to obtain the solution for the HE/EH, TE and TM modes, i.e. the full
% vector exact solution, in the three-layer structure. Unfortunately, while
% the published results are fine for TE and TM modes, the equations for the
% hybrid modes probably contains some typo and do not lead to the correct
% numerical solution. In 1997 ERDOGAN has approached the same problem once
% again and published the solution for the cladding modes (when light is
% guided by the cladding-surrounding interface) and suggested using the
% two-layer LP solution for the core-guided modes. His paper also contained
% a number of typos, and the errata was published in 2000. Finally, in 2005
% ZHANG and SHI published the full solution using the same approach as TSAO
% (1989). Unfortunately, this solution again seems to contain a mistake. 
% 
% The above relates to the eigen-value equation used to calculate the
% effective refractive index (or the propagation constant) of the guided
% mode. For full simulation, the light field (E and H vectors) should also
% be calculated. This is not implemented yet due to some problems with the
% references. BELANOV (1976) and MONERIE (1982) do not give the explicit
% solutions. ERDOGAN (1997, 2000) only gives it for the cladding-guided
% modes. ZHANG (2005) was not yet checked.
%
%% Simulation of n_eff in the tapered microfibre
% I suggest the following approach to simulation of the three-layer modes
% of the tapered microfibres. Initially, the fibre is weekly guiding and
% the Monerie solution is appropriate to calculate the effective refractive
% index in the core (will be proven below). As soon as the effective
% refractive index reaches the refractive index of the cladding, Erdogan
% solution is applied. In the region where the effective refractive index
% approaches that of the surrounding, this solution gets close to the
% two-layer solution --- classical HE/EH, TE and TM (proven below).
% modes), but there is no need to do so.
%
% Let's first calculate the exact two-layer modes to use them as a
% reference.
%
% We start by specifying the three materials in the order core, cladding,
% surrounging. As usual, the materials can be specified with refractive 
% indices or names known to |refrIndex| function. In the latter case,
% material dispersion is automatically taken into account. In the former
% case it is ignored.
materials = {1.455, 1.45, 'air'};
fibre.materials = materials;
fibre.coreCladdingRatio = 5/125;

% Specify the simulation task, first for the core
task = struct(...
    'nu', [0 1],... % first modal index
    'type', {{'hybrid', 'te', 'tm'}},... % mode types
    'maxmode', Inf,... % how many modes of each type and NU to calculate
    'lambda', 900,...
    'region', 'core');

% In the microfibre taper, the diameter varies from the diameter of the
% microwaist upto the diameter of the standard untapered fibre (125 um).
% All fibre dimensions in OFT are specified in micrometers, all wavelength
% values in nanometers.
argument = struct(...
    'type', 'dia',... % calculate mode curve n_eff vs. fibre diameter
    'harmonic', 1,... % required 
    'min', 0.01,... % minimum diameter
    'max', 125); % maximum diameter

addpath('..'); % path to OFT functions
infomode = false; % if true, OFT functions will output more information and 
% illustrate the simulation process with more pictures

%%
% Now calculate the mode curves:
if calculatemodes
modesTwoCore = buildModes(argument, fibre, task, infomode);
end
%% 
% Which mode was found?
t = modeDescription(modesTwoCore, false, false);
if iscell(t) 
    for i = 1:numel(t)-1, fprintf(' %s,', t{i}), end, fprintf(' %s\n', t{end});
else
    fprintf(' %s\n', t);
end
%%
% This is a so-called single-mode fibre
% similar to the well-known Fibercore SM800: in the core it guides only the
% fundamental mode HE_11 for all wavelengths above approximately 800 nm.   
%
% Display calculated curve
hTwoCore = showModes(modesTwoCore);
set(hTwoCore, 'Color', 'black')
%%
% Now repeat the same for the cladding modes. There is a lot of modes that
% can be guided in the cladding with 125 um diameter. To calculate all of
% them would take quite a bit of time, so we will calculate only the first 
% mode of each type (hybrid, TE and TM).
%
% We can use the same task specified above, just change it a little.
task.region = 'cladding';
task.maxmode = 1; % Calculate only the first of all families of modes, i.e.
% HE11, TE01 and TM01
if calculatemodes
    modesTwoClad = buildModes(argument, fibre, task, false);
end
hTwoClad = showModes(modesTwoClad);
set(hTwoClad, 'Color', 'black')
set(hTwoClad(2), 'LineStyle', '--') % Make TE curve dashed
%%
% Which modes have been calculated?
t = modeDescription(modesTwoClad, false, false);
fprintf('Calculated cladding two-layer modes: ');
if iscell(t) 
    for i = 1:numel(t)-1, fprintf(' %s,', t{i}), end, fprintf(' %s\n', t{end});
else
    fprintf(' %s\n', t);
end
%%
% The core modes are basically invisible now because 
% 
% $$n_\textrm{{\vphantom{lg}}core}-n_\textrm{clad}\ll n_\textrm{clad}-n_\textrm{surrounding}$$
% 
% We can resolve the curves by zooming vertically:
ylim([1.445 1.455])
%%
% The HE11 mode in the core does not reach the diameter of zero microns,
% due to numerical limitations. Theoretically, HE11 mode is guided at any
% fibre diameter, however small. 
%
% The TE_01 and TM_01 cladding modes are almost coinciding. We can
% zoom horizontally:
ylim auto
xlim([0 10])

%% Monerie modes
% Let's now calculate the Monerie modes in this fibre. In a separate demo
% on Monerie modes, I show that they should be traced from the core. We are
% interested in the LP01 mode corresponding to the HE11 mode as well as in
% LP11 corresponding to TE01 and TM01 modes. Our fibre core is single-mode,
% the LP11 mode is not supported. In order to trace it, we artificially
% increase the maximum considered diameter so that this mode is found and
% traced. Note: if argument.max is set to 180, there is a mode jump as
% explained in the Monerie tutorial
task = struct(...
    'nu', [0 1 2],... % first modal index
    'type', {{'monerie'}},... % mode types
    'maxmode', Inf,... % how many modes of each type and NU to calculate
    'lambda', 900,...
    'region', 'core');
argument.max = 200;
if calculatemodes
    modesMonerieCore = buildModes(argument, fibre, task, false);
end
hMonerieCore = showModes(modesMonerieCore);
set(hMonerieCore, 'Color', [1 0.5 0.5], 'LineWidth', 5)
set(hMonerieCore(2), 'LineStyle', '--')
uistack(hTwoCore, 'top')
uistack(hTwoClad, 'top')
xlim([0 125])
% Which modes have been found?
t = modeDescription(modesMonerieCore, false, false);
fprintf('Monerie (nu=[0 1 2]) modes found in the cladding: ');
if iscell(t) 
    for i = 1:numel(t)-1, fprintf(' %s,', t{i}), end, fprintf(' %s\n', t{end});
else
    fprintf(' %s\n', t);
end
%%
% As mentioned before, the Monerie solution is derived in the scalar
% approximation valid for small refractive index steps. Therefore we can
% expect it to coincide with the exact solution (hybrid, TE and TM modes)
% in the regions where the low refractive index surrounding does not play
% much role. This is the case at the large diameter because most of the
% field is still inside glass. At the same time, unlike the two-layer HE,
% TE and TM solutions, the Monerie solution gives a smooth transition
% between the core- to the cladding-guidance regions:
ylim ([1.449 1.451])
xlim([0 125])

%%
% For thin fibres, in the waist, the high refractive index step at the
% outer surface (between cladding and surrounding) makes the scalar Monerie
% solution for LP01 mode fully invalid, it strongly deviates from the HE11
% mode: 
ylim([1 1.45])
xlim([0 3])
%%
% The Monerie LP11 curve coincides, even for small diameters, with the TE01
% mode. This can be explained by the fact that a linear polarisation
% approximation is quite valid for the TE01 mode, which is 
% transversal for electric field. For TM01 mode, this approximation is
% not valid.
% 
% We have found that the Monerie solution is well suitable for
% the large diameter regions. We need a solution, which would work in the
% small diameter region and at the same time take into account the
% three-layer structure (so smoothly connect with the Monerie solution). 

%% Tsao modes
% As I mentioned before, the first known to me publication to treat this
% task was TSAO 1989. The eigen-value equation for the HE/EH modes is
% erroneous. But the TE and TM equations are fine.
%
% Let's calculate the TE01 and TM01 Tsao modes and see if they coincide
% with the two-layer T*01 solutions at small diameters and with the Monerie
% LP11 solution at the large diameter.
task = struct(...
    'nu', [0],... % first modal index
    'type', {{'tsaote', 'tsaotm'}},... % mode types
    'maxmode', 1,... 
    'lambda', 900,...
    'region', 'cladding');
argument.max = 125;
if calculatemodes
modesTsaoTClad = buildModes(argument, fibre, task, false);
end
% for i = 1:numel(modesTsaoTClad)
%     modesTsaoTClad(i) = addPointsToMode(modesTsaoTClad(i), [0 3]); 
%     modesTsaoTClad(i) = addPointsToMode(modesTsaoTClad(i), [0 3]); 
% end
hTsaoTClad = showModes(modesTsaoTClad);
set(hTsaoTClad, 'Color', 'green', 'LineWidth',3)
set(hTsaoTClad(1), 'LineStyle', '--')
uistack(hTwoCore, 'top')
uistack(hTwoClad, 'top')
ylim([1 1.45])
xlim([0 3])
%%
% As we see, at small diameters the Tsao TE and TM modes coincide well with
% the two-layer TE01 and TM01 modes. At larger diameter, the Tsao modes
% deviate from the two-layer solution and follow the Monerie LP11 mode, as
% expected.
ylim([1.44997 1.449975])
xlim([123 125])

%% Erdogan mode
% The three-layer fundamental mode (HE11) can be calculated using the
% solultion of ERDOGAN (1997, 2000).
task = struct(...
    'nu', 1,... % first modal index
    'type', {{'erdogan'}},... % mode types
    'maxmode', 1,... 
    'lambda', 900,...
    'region', 'cladding');
if calculatemodes
modesErdoganClad = buildModes(argument, fibre, task, false);
end
% modesErdoganClad = addPointsToMode(modesErdoganClad, [0 3]); 
% modesErdoganClad = addPointsToMode(modesErdoganClad, [0 3]); 
hErdoganClad = showModes(modesErdoganClad);
set(hErdoganClad, 'Color', 'cyan', 'LineWidth',3)
uistack(hTwoClad, 'top')
%%
% The found Erdogan mode nicely follows the two-layer HE11 mode in the
% small diameter region...
ylim auto
xlim([0 3])
%%
% ...and then deviates from it and "goes into the core", following the
% Monerie LP01 mode.
ylim ([1.4485 1.4505])
xlim([0 60])

%% Conclusions
% In this tutorial I have shown how to calculate the modal curves (n_eff
% vs. d) for three-layer fibres using OFT. As the example system, I've used
% a tapered optical microfibre. Calculation can be done using both
% two-layer and three-layer solutions.
%
% The two-layer solutions are exact in the regions where the third layer
% can be neglected and cannot simulate the transition region.
%
% The scalar approximation-based Monerie solution is valid for three-layer
% fibres when light 
% is confined inside glass so that it "does not see" the high
% refractive index step at the outer surface. This solution is therefore
% not valid for small diameters, at which a significant portion of light
% propagates in the evanescent field, but can be used to consider the
% transition region, where the two-layer solution does not provide enough
% information.
%
% For small diameter three-layer fibres, the Tsao and Erdogan
% solutions can be used. The Erdogan solution is only available for the
% cladding-guided modes (when n_eff < n_cladding). Therefore it seems
% reasonable to use the simpler Monerie solutions as long as the mode is
% core-guided (n_eff > n_cladding) and switch to the exact Erdogan and Tsao
% modes for n_eff < n_cladding.
%
%% Acknowledgements
% I am thankful to Timothy Lee and Peter Horak from ORC in Southampton,
% Ariane Stiebeiner from AG Rauschenbeutel in Vienna, and Fabian Bruse, our
% former diploma student, for their help during investigation of the
% available three-layer solutions. 
% 
%% References
% # Belanov et al., 1976: http://dx.doi.org/10.1070/QE1976v006n01ABEH010808
% # Erdogan, 1997: http://dx.doi.org/10.1364/JOSAA.14.001760<br/>
% # Erdogan, 2000: http://dx.doi.org/10.1364/JOSAA.17.002113
% # Monerie, 1982: http://dx.doi.org/10.1109/JQE.1982.1071586
% # Tsao et al., 1989: http://dx.doi.org/10.1364/JOSAA.6.000555
% # Zhang, Shi, 2005: http://dx.doi.org/10.1364/JOSAA.22.002516

save tutorial3ls.mat modesMonerieCore modesTsaoTClad modesTwoClad modesTwoCore modesErdoganClad