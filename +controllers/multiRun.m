function [results] = multiRun(parms)

% remove output of live graphics, does not make sense for multiple runs
parms.visualize = false;

workers = min(matlabpool('size'), parms.nRuns);

if (workers > 0)
    steps = floor(parms.nRuns / workers);
else
    steps = parms.nRuns;
    workers = 1;
end

if (workers*steps ~= parms.nRuns)
    tools.printLog(mfilename, 'WARNING: number of repetitions (%d) is not divisible by number of parallel workers (%d)!\n', parms.nRuns, workers);
    tools.printLog(mfilename, 'WARNING: reducing repetitions to %d!\n', workers*steps);
    parms.nRuns = workers*steps;
end

timings = zeros(workers,1);

if (exist('state', 'dir')==0)
    mkdir('.','state');
end
save('state/settings.mat','parms','steps','workers','-v7.3');

for k = 1:steps
    if exist('MDAPACK.stopcomputation','file')
        tools.printLog(mfilename, 'INFO: found "MDAPACK.stopcomputation" -  exiting repetitive loop and computing statistics\n', workers*steps);
        delete('MDAPACK.stopcomputation');
        % set repetitions to the actually performed repetitions
        parms.nRuns = workers*(k-1);
        % set steps to the actually taken steps
        steps = (k-1); %#ok<NASGU>
        
        save('state/settings.mat','steps','parms','-append');
        break;
    end
    
    try
        parfor i = 1:workers
            if (parms.progress) %#ok<PFBNS>
                tools.printLog(mfilename, 'repetition %3d of %3d: ', (k-1)*(workers)+i, parms.nRuns);
            end
            
            tmpres = [];
            retries = 10;
            finished = false;
            
            while (retries > 0 && ~finished)
                start = tic;
                try
                    tmpres = controllers.singleRun(parms);
                    finished = true;
                catch ex
                    if (strcmp(ex.identifier, 'MDAPACK:invalidRepresentation') && retries > 0)
                        % the ensemble was invalid during an assimilation run -
                        % we will retry the run
                        retries = retries - 1;
                        tools.printLog(mfilename, 'caught an MDAPACK:invalidRepresentation exception, retrying...\n');
                        continue;
                    else
                        % either we did not know the exception or the
                        % number of retries was exceeded -> propagate
                        rethrow(ex);
                    end
                end
                timings(i) = toc(start); %#ok<PFOUS>
            end
            
            res(i) = tmpres; %#ok<NASGU,PFOUS>
        end
        
        if ~exist('state', 'dir')
            mkdir('state');
        end
        
        % save the intermediate results
        save(sprintf('state/step_%05d.mat', k),'res','timings','k','-v7.3');
        save('state/settings.mat','k','-append','-v7.3');
        
        % free up the memory -- postprocessing is done later from the files
        clear('res');
        
    catch exception
        if (strcmp(exception.identifier, 'MATLAB:catenate:dimensionMismatch'))
        else
            % Unknown error. Just let it propagate.
            tools.printLog(mfilename, 'unknown error occured: %s ', exception.identifier);
            rethrow(exception);
        end
    end
end

results = controllers.postprocess();
end