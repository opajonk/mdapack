classdef DAPACKWrapper < models.Model
    %DAPACKWrapper DAPACK connection wrapper model
    %   The DAPACKWrapper model wraps around models implemented in DAPACK.
    
    properties
        config;
        evaluator;
        storage;
        tStep;
        decorrLength;
    end
    
    methods (Access = private)
        function H = computeH(this)
        end
        
        function R = computeR(this)
        end
        
        function L = computeL(this)
        end
        
        
        function [measurementAmout] = computeMeasurementAmount(this)
            measurements = this.config.getMeasurement();
            measurementAmout = 0;
            
            for i = 1:measurements.size
                var = this.storage.getDataVariance(measurements.get(i-1).getName(), this.tStep);
                if (var(1) > 0.0)
                    measurementAmout = measurementAmout + measurements.get(i-1).getLength();
                end
            end
        end
    end
    
    methods
        function this = DAPACKWrapper(configFileName, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            p.addRequired('configFileName', @(x) isa(x, 'char'));
            p.addRequired('pathToDAPACK', @(x) isa(x, 'char'));
            p.addRequired('nameServiceRef', @(x) isa(x, 'char'));
            p.addParamValue('decorrelationLength', 1.0, @(x)x>0);
            
            p.parse(configFileName, varargin{:});
            
            this.decorrLength = p.Results.decorrelationLength;
            
            tools.javalink.setdynpath(p.Results.pathToDAPACK);
            [this.config, this.evaluator, this.storage] = tools.javalink.dapack_link(p.Results.configFileName, p.Results.pathToDAPACK, p.Results.nameServiceRef);
            
            ensSize = this.config.getNumberOfExperiments();
            dvs = this.config.getDecisionVariables();
            
            % determine state vector length
            svLength = 0;
            for i = 1:dvs.getParameter().size
                svLength = svLength + dvs.getParameter().get(i-1).getLength();
            end
            
            for i = 1:dvs.getStateVariable().size
                svLength = svLength + dvs.getStateVariable().get(i-1).getLength();
            end
            
            % the measurements are appended to the statevector to enforce linearity of the measurement operator
            for i = 1:this.config.getMeasurement().size
                svLength = svLength + this.config.getMeasurement().get(i-1).getLength();
            end
            % reserve space for the ensemble
            this.ensemble = zeros(svLength, ensSize);
            
            % just allocate space for a "truth", left at 0 at all times
            this.truth = zeros(svLength, 1);
            
            % just make sure it has some value
            this.tStep = 1;
        end
        
        function step(this, t)
            ensSize = this.config.getNumberOfExperiments();
            dvs = this.config.getDecisionVariables();
            ms = this.config.getMeasurement();
            % copy matrix data to storage
            
            if (t > 1) % first timestep: we have nothing to update in the storage
                for p = 1:ensSize
                    % jump to start of "real" state vector
                    pointer = ms.size;
                    
                    for i = 1:dvs.getParameter().size
                        len = dvs.getParameter().get(i-1).getLength();
                        this.storage.setUpdatedDecisionVariable(dvs.getParameter().get(i-1).getName(), t-1, p, this.ensemble(pointer+1:pointer+len,p));
                        pointer = pointer + len;
                    end
                    
                    for i = 1:dvs.getStateVariable().size
                        len = dvs.getStateVariable().get(i-1).getLength();
                        this.storage.setUpdatedDecisionVariable(dvs.getStateVariable().get(i-1).getName(), t-1, p, this.ensemble(pointer+1:pointer+len,p));
                        pointer = pointer + len;
                    end
                    
                end
                clear pointer len i p;
            end
            
            
            % advance the ensemble
            this.evaluator.performComputation(t);
            
            % copy results from storage to matrix
            for p = 1:ensSize
                pointer = 0;
                
                % reset all measurements to 0.0
                this.ensemble(1:ms.size,p) = zeros(ms.size,1);
                
                % write measurements into statevector
                for i = 1:ms.size
                    len = ms.get(i-1).getLength();
                    if (len > 1); error('measurements of length > 1 are unsupported'); end
                    var = this.storage.getDataVariance(ms.get(i-1).getName(), t);
                    
                    if (var(1) > 0.0)
                        this.ensemble(pointer+1:pointer+len,p) = this.storage.getSimulatedData(ms.get(i-1).getName(), t, p);
                        pointer = pointer + len;
                    end
                end
                
                % jump to start of "real" state vector
                pointer = ms.size;
                
                for i = 1:dvs.getParameter().size
                    len = dvs.getParameter().get(i-1).getLength();
                    this.ensemble(pointer+1:pointer+len,p) = this.storage.getSimulatedDecisionVariable(dvs.getParameter().get(i-1).getName(), t, p);
                    pointer = pointer + len;
                end
                
                for i = 1:dvs.getStateVariable().size
                    len = dvs.getStateVariable().get(i-1).getLength();
                    this.ensemble(pointer+1:pointer+len,p) = this.storage.getSimulatedDecisionVariable(dvs.getStateVariable().get(i-1).getName(), t, p);
                    pointer = pointer + len;
                end
                
                
            end
            clear pointer len i p;
            
            % advance the truth: here, there is nothing to do as we already have some predefined history
            this.tStep = t;
        end
        
        function [err] = stepErr(this, t)
            assert(false, 'not implemented');
        end
        
        function [hasM] = hasMeasurement(this)
            hasM = (this.computeMeasurementAmount() > 0);
        end
        
        function [m] = measure(this)
            measurements = this.config.getMeasurement();
            length = this.computeMeasurementAmount();
            m = zeros(length, 1);
            
            pointer = 0;
            for i = 1:measurements.size
                len = measurements.get(i-1).getLength();
                if (len > 1); error('measurements of length > 1 are unsupported'); end
                var = this.storage.getDataVariance(measurements.get(i-1).getName(), this.tStep);
                
                if (var(1) > 0.0)
                    m(pointer+1:pointer+len,1) = this.storage.getDataMeasurement(measurements.get(i-1).getName(), this.tStep);
                    pointer = pointer + len;
                end
            end
            clear pointer len i;
        end
        
        function ret = ensMeasure(this)
            length = this.computeMeasurementAmount();
            ensSize = this.config.getNumberOfExperiments();
            measurements = this.config.getMeasurement();
            ret = zeros(length, this.ens_size);
            
            for p = 1:ensSize
                pointer = 0;
                for i = 1:measurements.size
                    len = measurements.get(i-1).getLength();
                    if (len > 1); error('measurements of length > 1 are unsupported'); end
                    var = this.storage.getDataVariance(measurements.get(i-1).getName(), this.tStep);
                    
                    if (var(1) > 0.0)
                        ret(pointer+1:pointer+len,p) = this.storage.getSimulatedData(measurements.get(i-1).getName(), this.tStep, p);
                        pointer = pointer + len;
                    end
                end
                clear pointer len i;
            end
        end
        
        function [H] = measureOp(this)
            H = zeros(size(this.ensemble, 1), this.computeMeasurementAmount());
            ms = this.config.getMeasurement();
            k = 1;
            
            for i = 1:ms.size
                len = ms.get(i-1).getLength();
                if (len > 1); error('measurements of length > 1 are unsupported'); end
                var = this.storage.getDataVariance(ms.get(i-1).getName(), this.tStep);
                
                if (var(1) > 0.0)
                    H(k,k) = 1.0;
                    k = k + 1;
                end
            end
            
            H = H';
        end
        
        function [merr] = measureErr(this)
            length = this.computeMeasurementAmount();
            merr = (chol(this.measureCov()) * randn(length,1));
        end
        
        function [cov] = measureCov(this)
            measurements = this.config.getMeasurement();
            length = this.computeMeasurementAmount();
            cov = zeros(length);
            
            pointer = 0;
            for i = 1:measurements.size
                len = measurements.get(i-1).getLength();
                if (len > 1); error('measurements of length > 1 are unsupported'); end
                
                var = this.storage.getDataVariance(measurements.get(i-1).getName(), this.tStep);
                
                if (var(1) > 0.0)
                    cov(pointer+1:pointer+len,pointer+1:pointer+len) = diag(var(1));
                    pointer = pointer + len;
                end
            end
        end
        
        function [L] = distanceMatrix(this)
            ordering = javaArray('java.lang.String[]',2);
            
            ms = this.config.getMeasurement();
            ordering1 = javaArray('java.lang.String', this.computeMeasurementAmount());
            pos = 1;
            for i = 1:ms.size
                var = this.storage.getDataVariance(ms.get(i-1).getName(), this.tStep);
                if (var(1) > 0.0)
                    ordering1(pos) = ms.get(i-1).getName();
                    pos = pos + 1;
                end
            end
            
            dvs = this.config.getDecisionVariables();
            ordering2 = javaArray('java.lang.String', dvs.getParameter().size() + dvs.getStateVariable().size());
            pos = 1;
            for i = 1:dvs.getParameter().size
                ordering2(pos) = dvs.getParameter().get(i-1).getName();
                pos = pos + 1;
            end
            
            for i = 1:dvs.getStateVariable().size
                ordering2(pos) = dvs.getStateVariable().get(i-1).getName();
                pos = pos + 1;
            end
            
            ordering(1) = ordering1;
            ordering(2) = ordering2;
            
            L = this.evaluator.getMethodParameter(java.lang.String('ensembleSensorDistances'), this.tStep, ordering);
        end
        
        function [dl] = decorrelationLength(this)
            dl = this.decorrLength;
        end
        
        function [tCount] = timeStepCount(this)
            tCount = this.config.getTimeSteps();
        end
        
        function plot(this, varargin)
            % ??
        end
        
        function [this] = plotMeasurementPos(this)
            % ??
        end
        
        function [this] = plotEnsemble(this)
            % ??
        end
        
        %% special functions for localization with wavelets
        
        function [gridSize] = getGridSize(this)
            p = this.config.getDecisionVariables().getParameter();
            gridSize = this.evaluator.getMethodParameter('gridSize', 1, p.get(1).getName());
            
            assert(length(gridSize) == 3);
        end
        
        % return the wavelet transformed ensemble and some helper variables
        % WARNING: HACK HACK HACK HACK, not a nice solution!
        function [wX, S, decompositionLength] = wavedec(this, X, waveletLevel, wavelet)
            [gridSize, ensSize, cellCount, activeCellCount, gridVariableCount, gridMapping] = this.computeHelperVariables();
            dvs = this.config.getDecisionVariables();
            ms = this.config.getMeasurement();
            totalVariableCount = dvs.getParameter().size()  + dvs.getStateVariable().size();
            scalarVariableCount = totalVariableCount - gridVariableCount;
            
            % compute length of the wavelet transform of one z-level of a single grid variable
            tmp = this.storage.getSimulatedDecisionVariable(dvs.getParameter().get(0).getName(), this.tStep, 1);
            tmp = reshape(tmp, gridSize');
            [C,S] = wavedec2(tmp(:,:,1),waveletLevel,wavelet); 
            decompositionLength = length(C);
            clear tmp C;
            
            % reserve memory for wavelet transformed ensemble and copy measurements
            wX = zeros([decompositionLength * gridSize(3) * gridVariableCount + scalarVariableCount + ms.size(), ensSize]);
            wX(1:ms.size(),:) = X(1:ms.size(),:);
            
            for p = 1:ensSize
                % jump to start of "real" state vector (behind measurements)
                pointerX = ms.size();
                pointerWX = ms.size();
                
                % Store wavelet transformed grid parameters in the state vector. Those which originate from
                % a grid are reshaped to the 3D grid and then wavelet transformed. The other ones, e.g. scalars,
                % are simply copied. This is done for parameters and state variables alike. The latter ones are mapped to a
                % full grid to be able to compute the wavelet transform consistently.
                for i = 1:dvs.getParameter().size
                    len = dvs.getParameter().get(i-1).getLength();
                    
                    if (len == cellCount) % we have a grid-based variable which we have to wavelet transform
                        dataTmp = X(pointerX+1:pointerX+cellCount, p);
                        dataTmp = reshape(dataTmp, gridSize');
                        
                        for z = 1:gridSize(3) % for each z-level we do a separate 2D wavelet transform
                            wX(pointerWX+1:pointerWX+decompositionLength,p) = wavedec2(dataTmp(:,:,z),waveletLevel,wavelet);
                            pointerWX = pointerWX + decompositionLength;
                        end
                    else % other variable (possibly scalar), so just copy it
                        wX(pointerWX+1:pointerWX+len,p) = X(pointerX+1:pointerX+len,p);
                        pointerWX = pointerWX + len;
                    end
                    
                    pointerX = pointerX + len;
                end
                
                
                for i = 1:dvs.getStateVariable().size
                    len = dvs.getStateVariable().get(i-1).getLength();
                    if (len == activeCellCount) % we have a grid-based variable which we have to wavelet transform
                        dataTmp = X(pointerX+1:pointerX+activeCellCount, p);
                        dataTmp = gridMapping * dataTmp; % map from active to full grid
                        dataTmp = reshape(dataTmp, gridSize');
                        
                        for z = 1:gridSize(3) % for each z-level we do a separate 2D wavelet transform
                            wX(pointerWX+1:pointerWX+decompositionLength,p) = wavedec2(dataTmp(:,:,z),waveletLevel,wavelet);
                            pointerWX = pointerWX + decompositionLength;
                        end
                    else % other variable (possibly scalar), so just copy it
                        wX(pointerWX+1:pointerWX+len,p) = X(pointerX+1:pointerX+len,p);
                        pointerWX = pointerWX + len; 
                    end
                    
                    pointerX = pointerX + len;
                end
            end
        end
        
        % reconstruct the ensemble from wavelet transform obtained by wavedec and store it
        function [X] = waverec(this, wavelet, wX, S, decompositionLength)
            [gridSize, ensSize, cellCount, activeCellCount, ~, gridMapping] = this.computeHelperVariables();
            dvs = this.config.getDecisionVariables();
            ms = this.config.getMeasurement();
            
            % reserve memory and copy measurements
            X = zeros(size(this.ensemble));
            X(1:ms.size(),:) = wX(1:ms.size(),:);
            
            % perform a 2D wavelet reconstruction for each z-level of each grid variable and ensemble member, copy other variables
            for p = 1:ensSize
                % start after the measurements
                pointerX = ms.size();
                pointerWX = ms.size();
                
                for i = 1:dvs.getParameter().size
                    len = dvs.getParameter().get(i-1).getLength();
                    
                    if (len == cellCount) % we have a grid based parameter
                        dataTmp = zeros(gridSize');
                        
                        for z = 1:gridSize(3) % perform a 2D wavelet reconstruction for each z-level
                            dataTmp(:,:,z) = waverec2(wX(pointerWX+1:pointerWX+decompositionLength,p),S,wavelet);
                            pointerWX = pointerWX + decompositionLength;
                        end
                        
                        % reshape 3D data to column vector
                        dataTmp = reshape(dataTmp,[],1);
                        % store to ensemble matrix
                        X(pointerX+1:pointerX+cellCount,p) = dataTmp;
                    else % no grid based variable, so just copy it
                        X(pointerX+1:pointerX+len,p) = wX(pointerWX+1:pointerWX+len,p);
                        pointerWX = pointerWX + len;
                    end
                    
                    pointerX = pointerX + len;
                end
                
                for i = 1:dvs.getStateVariable().size
                    len = dvs.getStateVariable().get(i-1).getLength();
                    
                    if (len == activeCellCount)
                        dataTmp = zeros(gridSize');
                        
                        for z = 1:gridSize(3) % perform a 2D wavelet reconstruction for each z-level
                            dataTmp(:,:,z) = waverec2(wX(pointerWX+1:pointerWX+decompositionLength,p),S,wavelet);
                            pointerWX = pointerWX + decompositionLength;
                        end
                        
                        % reshape to column vector
                        dataTmp = reshape(dataTmp,[],1);
                        % map to active grid
                        dataTmp = gridMapping' * dataTmp;
                        % store to ensemble matrix
                        X(pointerX+1:pointerX+activeCellCount,p) = dataTmp;
                    else % no grid based variable, so just copy it
                        X(pointerX+1:pointerX+len,p) = wX(pointerWX+1:pointerWX+len,p);
                        pointerWX = pointerWX + len;
                    end
                    
                    pointerX = pointerX + len;
                end
            end
        end
        
        function [gridSize, ensSize, cellCount, activeCellCount, gridVariableCount, gridMapping] = computeHelperVariables(this)
            % determine a buch of helper variables and data from the model
            gridSize = this.evaluator.getMethodParameter('gridSize', 1, this.config.getDecisionVariables().getParameter().get(1).getName());
            assert(length(gridSize) == 3);
            activeGridMap = this.evaluator.getMethodParameter('activeGridMap', 1, this.config.getDecisionVariables().getStateVariable().get(1).getName());
            dvs = this.config.getDecisionVariables();
            
            ensSize = this.config.getNumberOfExperiments();
            cellCount = length(activeGridMap);
            activeCellCount = nnz(activeGridMap);
            
            % compute the mapping operator between active grid and full grid
            a = 1:activeCellCount;
            b = zeros(1,activeCellCount);
            pointer = 1;
            
            for i = 1:cellCount
                if(activeGridMap(i) == true)
                    b(pointer) = i;
                    pointer = pointer + 1;
                end
            end
            
            % identify the amount of parameters/state variables which actually represent grid data
            % (as opposed to scalar meta parameters for relperms, for example)
            % Only these will be wavelet-transformed, obviously!
            gridVariableCount = 0;
            for i = 1:dvs.getParameter().size
                if (dvs.getParameter().get(i-1).getLength() == cellCount)
                    gridVariableCount = gridVariableCount + 1;
                end
            end
            for i = 1:dvs.getStateVariable().size
                if (dvs.getStateVariable().get(i-1).getLength() == activeCellCount)
                    gridVariableCount = gridVariableCount + 1;
                end
            end
            
            % right multiply this operator with data from the active grid to map it to the full grid
            gridMapping = sparse(b,a,ones(1,activeCellCount),cellCount,activeCellCount,activeCellCount);
        end
    end
end

