function results = batchcall(parameters, func, resultFilePattern, varargin)
%BATCHCALL Easily call a MATLAB function with multiple parameter sets

[experiments, varnames] = determineExperiments(parameters);

for experimentNo = 1:size(experiments, 1)
    for j = 1:length(varnames)
        paramSet(experimentNo).(varnames{j}) = experiments{experimentNo, j}; %#ok<AGROW>
    end
end

for experimentNo = 1:size(experiments, 1)
    parameters = paramSet(experimentNo);
    
    if (parameters.progress)
        tools.printLog(mfilename, 'performing experiment %d of %d\n', experimentNo, size(experiments, 1));
    end
    
    machine = computer(); %#ok<NASGU>
    
    startTimestamp = datestr(now()); %#ok<NASGU>
    results = func(parameters, varargin{:});
    endTimestamp = datestr(now()); %#ok<NASGU>
    
    if exist('resultFilePattern', 'var')
        if ~exist('results', 'dir')
            mkdir('results');
        end
        filename = createFileName(resultFilePattern, parameters, experimentNo);
        save(filename, '-v7.3', 'parameters', 'results', 'experimentNo', 'startTimestamp', 'endTimestamp', 'machine');
    end
    
    if (parameters.progress)
        tools.printLog(mfilename, 'end of experiment %d\n', experimentNo);
    end
end
end % function batchcall

function filename = createFileName(resultFilePattern, paramSet, experimentNo)
names = fieldnames(paramSet);
for j = 1:length(names)
    toReplace = strcat('\{', names(j), '\}');
    % TODO possibly we need other conversions here
    if (isnumeric(paramSet.(names{j})))
        replacement = num2str(paramSet.(names{j}));
    elseif (islogical(paramSet.(names{j})))
        if paramSet.(names{j})
            replacement = 'true';
        else
            replacement = 'false';
        end
    else
        replacement = char(paramSet.(names{j}));
    end
    resultFilePattern = regexprep(resultFilePattern, toReplace, replacement);
end

filename = sprintf(resultFilePattern, experimentNo);
end

function [ experiments, names ] = determineExperiments( parameters )
%determineExperiments Creates a matrix of all parameter combinations

% compute how many experiments we have to run by taking the product of all
% parameter lengths
expCount = 1;
for x = fieldnames(parameters)'
    expCount = expCount * length(parameters.(x{1}));
end

% create the experiments cell array containing all possible parameter
% combinations
experiments = cell(expCount, length(fieldnames(parameters)));
stride = expCount;
names = fieldnames(parameters);
for j = 1:length(names)
    varlen = length(parameters.(names{j}));
    stride = stride / varlen;
    values = parameters.(names{j});
    for i = 1:expCount
        idx = mod(floor((i-1)/stride), varlen)+1;
        assert(~isstruct(values{idx}), 'no structured values allowed');
        experiments{i,j} = values{idx};
    end
end
end

