function result = styleVsModetype(modetype, nu)
% Creates line-style understood by SET according to two-layer mode type

hybrid = strcmpi(modetype , 'HYBRID');
te = strcmpi(modetype , 'TE');
tm = strcmpi(modetype , 'TM');

if hybrid == 1
    if nu == 1
        result = ':'; 
    elseif nu == 2
        result = '-';
    else
        result = '--';
    end;
elseif te == 1
    result = '-.';
elseif tm == 1
    result = '-.';
else
    throw(MException('InvalidModetype', '%s: Invalid mode type (%s)\n', upper(mfilename), upper(modetype)));
end;