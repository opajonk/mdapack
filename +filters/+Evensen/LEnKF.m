classdef LEnKF < filters.Filter
    %LENKF Ensemble Kalman filter (EnKF) using perturbed observations and local updating
    
    properties (Access = private)
        opts;
        
        measurementEnsembleRng; % the PRNG used to generate the measurement ensemble
    end
    
    methods
        function this = LEnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('sampling', 'default', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            p.addOptional('samplingFactor', 1.0, @(x)x>=1.0 );
            p.addOptional('solve', 'eigenDecomp', ...
                @(x)any(strcmpi(x,{'eigenDecomp','MATLAB'})));
            p.addOptional('sampleBias', 'remove', ...
                @(x)any(strcmpi(x,{'keep','remove'})));
            p.addOptional('measurementCov', 'exact', ...
                @(x)any(strcmpi(x,{'exact','low-rank'})));
            
            %% option parsing
            p.parse(varargin{:});
            this.opts = p.Results;
            
            % extract the PRNG from the options - it is not an option, so
            % remove it from there
            this.measurementEnsembleRng = this.opts.measurementEnsembleRng;
            this.opts = rmfield(this.opts, 'measurementEnsembleRng');
        end
        
        function this = update(this, model, representation)
            % preparations
            [A sample] = representation.getEnsemble();
            M = size(A, 1);
            N = size(A, 2);
            L = model.distanceMatrix();
            dL = 2*model.decorrelationLength();
            L = arrayfun(@(x) (x <= dL) , L);
                        
            % get noisy measurement
            d = model.measure();
            
            % sample from measurement error
            if strcmp(this.opts.sampling, 'advanced') && this.opts.samplingFactor > 1.0
                E = model.measureErr(this.opts.samplingFactor*N, this.measurementEnsembleRng);                
                E = tools.sampling.condenseSample(E, N);
            else
                E = model.measureErr(N, this.measurementEnsembleRng); 
            end
            
            if strcmp(this.opts.sampleBias, 'remove')
                % OPTIONAL: Remove sampling bias
                E = tools.sampling.removeBias(E);
            end
            
            % create measurement ensemble
            D = bsxfun(@plus, E, d);
            
            % compute ensemble representation of measurement error
            % covariance matrix
%             Re = (1/(N-1)) * (E * E');
            
            % obtain measurements from the ensemble
            h = model.measureOp();
            HA = h(representation);
            
            % compute ensemble of innovation vectors
            Dprime = D - HA;
            
            % compute measurement mean
            meanHA = mean(HA,2);
            
            % compute ensemble mean
            meanA = mean(A, 2);
            
            % substract effect mean from effects on ensemble
            Aprime  =  A - repmat(meanA,  1, size(A,2));
            HAprime = HA - repmat(meanHA, 1, size(HA,2));
            
            if strcmp(this.opts.measurementCov, 'exact')
                % exact representation of the measurement error covariance
                % matrix
                C = HAprime*HAprime' + (N-1)*model.measureCov();
            elseif strcmp(this.opts.measurementCov, 'low-rank')
                % low-rank approximation to the measurement error
                % covariance (see Evensen2004, eq. (12) & (14) and
                % Evensen2009a for a discussion)
                C = HAprime*HAprime' + E*E';
            end
            
            % Now we have everything ready to compute the real update
            
            if strcmp(this.opts.solve, 'MATLAB')
                % Variant 1: Solve using MATLAB backslash operator
                
                for i=1:M
                    % Setup "chooser matrix" for localisation. This maps
                    % the space of all measurements to the space of
                    % "measurements relevant to grid point i".
                    Y = diag(L(i,:));
                    deleteRows = [];
                    for j = 1:length(Y)
                        if sum(Y(j,:)) == 0
                            deleteRows = [deleteRows j]; %#ok<AGROW>
                        end
                    end
                    Y(deleteRows, :) = [];
                    
                    % it is possible that no measurement is relevant for
                    % grid point i --> no update!
                    if (~isempty(Y))
                        A(i,:) = A(i,:) + Aprime(i,:) * (HAprime' * Y') * ((Y*C*Y') \ Y*Dprime);
                    end
                end
                
            elseif strcmp(this.opts.solve, 'eigenDecomp')
                % Variant 2: Solve using (limited) eigenvalue decomposition
                % like described in Evensen2003, eqn. (55,56).
                % Important: different types of measurements must be
                % assimilated sequentially, or scaled.
                                
                for i=1:M
                    % Setup "chooser matrix" for localisation. This maps
                    % the space of all measurements to the space of
                    % "measurements relevant to grid point i".
                    Y = diag(L(i,:));
                    deleteRows = [];
                    for j = 1:length(Y)
                        if sum(Y(j,:)) == 0
                            deleteRows = [deleteRows j]; %#ok<AGROW>
                        end
                    end
                    Y(deleteRows, :) = [];
                    
                    % it is possible that no measurement is relevant for
                    % grid point i --> no update!
                    if (~isempty(Y))
                        [V, D] = eig(Y*C*Y');
                        % [V, D] = eigs(HAprime*HAprime' + E*E',[], min(N,length(d)));
                        % eigs is a little bit slower than eig...
                        weights = V(1:min(N,size(Y,1)),1:min(N,size(Y,1))) * diag(1./diag(D(1:min(N,size(Y,1)),1:min(N,size(Y,1))))) * V(1:min(N,size(Y,1)),1:min(N,size(Y,1)))';
                        A(i,:) = A(i,:) + Aprime(i,:) * (HAprime' * Y') * weights * (Y*Dprime);
                    end
                end
            end
            
            representation.setEnsemble(A, sample);
        end
        
        function stat = getStatistics(this)
            stat = struct([]);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

