function [pmpsNew, harmonics] = setHarmonic(pmps)
% Sets the HARMONIC parameter of a phase-matching point structure.
% Also sets HARMONIC of both modes within this structure.

% Copyright: (cc-by) Subwavelength-diameter fibres team @ Uni Bonn, 2010 
% http://agmeschede.iap.uni-bonn.de, kotya.karapetyan@gmail.com

for i = 1:length(pmps)
    pmp = pmps(i);
    if strcmpi(pmp.argtype, 'WVL')
        harmonic = max(pmp.mode1.harmonic, pmp.mode2.harmonic);
        if abs(pmp.mode1.harmonic - pmp.mode2.harmonic) < 1e-6 % harmoincs are equal
            error('Cannot determine the harmonic value now\n');
        end;
    elseif strcmpi(pmp.argtype, 'DIA')
        if abs(mod(pmp.mode1.par, pmp.mode2.par)) < 1e-6
            harmonic = round(pmp.mode1.par / pmp.mode2.par);
            pmp.mode1.harmonic = 1;
            pmp.mode2.harmonic = harmonic;
        elseif abs(mod(pmp.mode2.par, pmp.mode1.par)) < 1e-6
            pmp.mode2.harmonic = 1;
            harmonic = round(pmp.mode2.par / pmp.mode1.par);
            pmp.mode1.harmonic = harmonic;
        else
            error('Cannot determine the harmonic');
        end
    else
        error('Invalid argtype');
    end
    
    if pmp.mode1.harmonic > pmp.mode2.harmonic
        % shift modes so that first mode is fundamental
        modeTemp = pmp.mode2;
        pmp.mode2 = pmp.mode1;
        pmp.mode1 = modeTemp;
        neffTemp = pmp.neff1;
        pmp.neff1 = pmp.neff2;
        pmp.neff2 = neffTemp;
        clear modeTemp;
        clear neffTemp;
    end;
    
    pmp.harmonic = harmonic; % set harmonic in pmp itself too
    harmonics(i) = harmonic; %#ok<AGROW>
    pmpsNew(i) = pmp; %#ok<AGROW>
end;