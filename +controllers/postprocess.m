function [results parms] = postprocess(deleteFiles) %#ok<*STOUT>

if (nargin < 1)
    deleteFiles = false;
else
    assert(islogical(deleteFiles))
end

load('state/settings.mat');

%% postprocessing results
if (parms.progress)
    tools.printLog(mfilename, 'postprocessing results... \n');
end

load(sprintf('state/step_%05d.mat', 1),'res','timings');

% initialize output structure
names = fieldnames(res(1)); %#ok<NODEF>
results.nRuns = parms.nRuns;
results.timings = [];

if (parms.statistics && parms.nRuns > 1)
    for j = 1:length(names)
        results.(names{j}) = [];
        results.(names{j}).min = [];
        results.(names{j}).max = [];
        results.(names{j}).mean = [];
        results.(names{j}).var = [];
    end
else
    for j = 1:length(names)
        results.(names{j}) = [];
    end
end

for k = 1:steps
    if (parms.progress)
        tools.printLog(mfilename, 'postprocessing step %d... \n', k);
        tic;
    end
    
    load(sprintf('state/step_%05d.mat', k),'res','timings');
    results.timings = [results.timings; timings];
    
    if (parms.statistics && parms.nRuns > 1)        
        tic;
        if (parms.progress)
            tools.printLog(mfilename, 'updating statistics with %d runs... \n', workers);
        end
        
        for j = 1:length(names)
            for i = 1:workers
                t = (k-1)*workers + i;
                
                if ( t <= 1)
                    results.(names{j}).min = res(i).(names{j});
                    results.(names{j}).max = res(i).(names{j});
                    
                else
                    results.(names{j}).min = min([results.(names{j}).min res(i).(names{j})], [], ndims(results.(names{j})));
                    results.(names{j}).max = max([results.(names{j}).max res(i).(names{j})], [], ndims(results.(names{j})));
                end
                
                % recursively update the mean, see
                % [WeissteinSampleVarianceComputation]
                if ( t <= 1)
                    results.(names{j}).mean = res(i).(names{j});
                else
                    mold = results.(names{j}).mean;
                    results.(names{j}).mean = results.(names{j}).mean + (res(i).(names{j}) - results.(names{j}).mean)./(t+1);
                end
                
                % recursively update the variance, see
                % [WeissteinSampleVarianceComputation]
                if ( t <= 1)
                    results.(names{j}).var = zeros(size(res(i).(names{j})));
                else
                    results.(names{j}).var = (1-(1/t)) .* results.(names{j}).var + (t+1) .*(results.(names{j}).mean - mold).^2;
%                     results.(names{j}).var = ((t-1)/t) .* results.(names{j}).var + (1/(t-1)) .*(res(i).(names{j}) - results.(names{j}).mean).^2;
                end
            end
        end
        
        if (parms.progress)
            tools.printLog(mfilename, 'done\n ');
            toc;
        end
    else
        if (parms.progress)
            tools.printLog(mfilename, 'concatenating results from %d runs... \n', workers);
        end
        
        for j = 1:length(names)
            catdim = ndims(res(1).(names{j}))+1;
            catcommand = ['cat(',int2str(catdim),','];
            
            % if these are not the first results, we have to prepend the
            % already computed ones
            if exist('results','var') && isstruct(results) && isfield(results, names{j})
                catcommand = [catcommand, 'results.', names{j} ,',']; %#ok<AGROW>
            end
            
            for i = 1:workers
                if (i == workers)
                    catcommand = [catcommand, 'res(',int2str(i),').(names{',int2str(j),'}));']; %#ok<AGROW>
                else
                    catcommand = [catcommand, 'res(',int2str(i),').(names{',int2str(j),'}), ']; %#ok<AGROW>
                end
            end
            
            results.(names{j}) = eval(catcommand);
            
            % free memory
            for i = 1:workers
                res(i).(names{j}) = []; %#ok<AGROW>
            end
        end
        if (parms.progress)
            tools.printLog(mfilename, 'done\n ');
            toc;
        end
    end
end

% optionally remove the state files
if (deleteFiles)
    delete('state/settings.mat');
    for k = 1:steps
        delete(sprintf('state/step_%05d.mat', k));
    end
end
end

