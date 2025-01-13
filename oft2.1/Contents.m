%% Optical Fibre Toolbox
% Optical Fibre Toolbox (OFT) provides functions for fast automatic
% calculation of guided modes in simple optical fibres. Developed with
% tapered microfibres (aka nanofibres) in mind. Exact solutions for weak
% and strong guidance cases are provided. Material dispersion is taken into 
% account.
%
% For authors, revisions see readme.txt.
% For EXAMPLES of use, see tutorial3ls.m and oftDemo.m in the DEMO directory.
%
% Main functions
%   buildModes          - Calculates the modal curves (n_eff vs d or n_eff vs lambda)
%   fibreForPM          - Returns possible phase-matching points for 2nd or 3rd harmonic generation
%   modeField           - Returns the field (E, H) value(s) for the two-layer modes (HE/EH, TE, TM)
%   neff                - Returns n_eff for particular mode and argument value
%
% Plotting functions
%   showModes           - Displays modal curves for the MODES array.
%   showPMPoints        - Marks phase-matching points (mode curve intersections) on the current plot 
%   displayField2       - Creates plots of the mode field for two-layer structure
%
% Additional functions
%   addPointsToMode     - Adds points to the mode curve, trying to double their density 
%   fieldPower          - Power in the mode field.
%   modeDescription     - Returns text description of the mode(s).
%   modeVsV             - Returns the same mode with V-parameter as the argument
%   normalizePower      - Returns field structure normalised to the given power
%   overlapIntegral     - Overlap integral for two modes at a phase-matching point 
%   phaseMatchingPoints - Finds all phase-matching points for the given modes
%   pmpDescription      - Creates a text description of a phase-matching point
%   refrIndex           - Returns refractive index of some known materials at a given wavelength
% 
% Technical functions 
%   besselid            - Derivative of the modified Bessel function of the first kind
%   besseljd            - Derivative of the Bessel function of the first kind
%   besselkd            - Derivative of the modified Bessel function of the second kind
%   besselyd            - Derivative of the Bessel function of the second kind
%   checkTask           - Verifies task correctness for the given fibre structure
%   colourVsLambda      - Returns a colour as corresponding to the given wavelength. 
%   correctModeOrder    - Checks that NU is allowed for such MODETYPE.
%   eve2LS              - Eigen-value equations for vector two-layer structure modes
%   eve3LS              - Eigen-value equations for three-layer structure modes
%   fibreMode           - General eigen-value equation function, calls eve2LS or eve3LS
%   fieldGrid           - Two-layer grid uniform in R and PHI, for displayField2
%   findModes           - Returns starting points n_eff(arg) for all modes found according to task
%   findRoot            - Works on top of fzero, tries to make it work in large number of cases.
%   modeCheckValue      - Check-sum for the modes array.
%   partype             - Returns the parameter type ('wvl' or 'dia')
%   phaseMatch          - Tries to find a phase-match point on a between the two given modes
%   pmCheckValue        - Copyright: (cc-by) Subwavelength-diameter fibres team @ Uni Bonn, 2010 
%   pmSameModes         - Returns TRUE if PMP1 and PMP2 are the intersections of the same two modes
%   setHarmonic         - Sets the HARMONIC parameter of a phase-matching point structure.
%   styleVsModetype     - Creates line-style understood by SET according to two-layer mode type
%   traceMode           - Calculates the mode curve (n_eff vs d or lambda) starting from one point
%   vParameter          - Returns the V-parameter for two-layer structures 
%
% Copyright (cc-by) K. Karapetyan et al., AG Meschede, Uni Bonn, 2008--2011
% kotya.karapetyan@gmail.com, http://agmeschede.iap.uni-bonn.de
% For a full list of contributors please see authors.txt. 
