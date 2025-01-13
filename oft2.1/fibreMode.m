function result = fibreMode(d, neff, lambda, fibreSpec, modeTask)
% General eigen-value equation function, calls eve2LS or eve3LS

switch lower(modeTask.modetype)
    case {'hybrid', 'te', 'tm'}
        result = eve2LS(d, neff, lambda, fibreSpec, modeTask);
    case {'monerie', 'tsaote', 'tsaotm', 'erdogan', 'zhang'}
        result = eve3LS(d, neff, lambda, fibreSpec, modeTask);
    otherwise
        error('Wrong mode type');
end
