function result = pmSameModes(pmp1, pmp2) 
% Returns TRUE if PMP1 and PMP2 are the intersections of the same two modes
% (with different parameters)

mode11 = modeDescription(pmp1.mode1, false);
mode12 = modeDescription(pmp1.mode2, false);
mode21 = modeDescription(pmp2.mode1, false);
mode22 = modeDescription(pmp2.mode2, false);

result = (strcmpi(mode11, mode21) && strcmpi(mode12, mode22)) ||...
    (strcmpi(mode11, mode22) && strcmpi(mode12, mode21));