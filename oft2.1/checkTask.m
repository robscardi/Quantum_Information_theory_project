function checkTask(task, fibre)
% Verifies task correctness for the given fibre structure

if iscell(task.type)
    numberOfModeTypes = length(task.type);
else
    numberOfModeTypes = 1;
end

if numberOfModeTypes > 1
    for i = 1:numel(numberOfModeTypes)
        singleModeTypeTask = task;
        singleModeTypeTask.type = task.type{i};
        checkTask(singleModeTypeTask, fibre);
    end
else    
    threeLayerMode = strcmpi(task.type, 'monerie') || strcmpi(task.type, 'erdogan') ...
        || strcmpi(task.type, 'tsaote') || strcmpi(task.type, 'tsaotm') ...
        || strcmpi(task.type, 'zhang');
    twoLayerMode = strcmpi(task.type, 'hybrid') || strcmpi(task.type, 'te') ...
        || strcmpi(task.type, 'tm');
    
    if ~threeLayerMode && ~twoLayerMode
        error('Wrong mode type');
    end
    
    if threeLayerMode
        assert(numel(fibre.materials) == 3, 'Three layer mode and ~3 layers');
    end
    
    if numel(fibre.materials) == 3
        assert(strcmpi(task.region, 'core') || strcmpi(task.region, 'cladding'),...
            'Core or cladding region expected');
    end
end