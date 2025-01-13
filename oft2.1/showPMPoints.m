function showPMPoints(PMPOINTS)
% Marks phase-matching points (mode curve intersections) on the current plot 

if isempty(PMPOINTS)
    return;
end;

for i=1:length(PMPOINTS)
    assert((PMPOINTS(i).quality == 1) || (PMPOINTS(i).quality == -1));
    if PMPOINTS(i).quality == 1
        colour = 'k';
    else
        colour = 'r';
    end;
    
    plot(PMPOINTS(i).arg, PMPOINTS(i).neff1, 'Marker', 'o', 'Color', colour, 'MarkerSize', 10);
end;
