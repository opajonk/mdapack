classdef EnKF < filters.Filter
    %ENKF Ensemble Kalman filter (EnKF) using perturbed observations
    
    properties (Access = private)
        opts;
        
        measurementEnsembleRng; % the PRNG used to generate the measurement ensemble
        currentRho;
    end
    
    methods
        function this = EnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('sampling', 'default', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            p.addOptional('samplingFactor', 1.0, @(x)x>=1.0 );
            p.addOptional('samplingOrder', 2, @(x) any(find([0,1,2] == x)));
            p.addOptional('solving', 'eigen', ...
                @(x)any(strcmpi(x,{'eigen','MATLAB'})));
            p.addOptional('sampleBias', 'remove', ...
                @(x)any(strcmpi(x,{'keep','remove'})));
            p.addOptional('measurementCov', 'exact', ...
                @(x)any(strcmpi(x,{'exact','low-rank'})));
            p.addOptional('inflation', 'off', ...
                @(x)any(strcmpi(x,{'off','fixed','adaptive'})));
            p.addOptional('gamma', 1.0, @(x)x>=1.0 );
            p.addOptional('fixedRho', 1.0, @(x)x>=1.0 );
            p.addOptional('measurementEnsembleRng', RandStream('mlfg6331_64'), @(x) isa(x,'RandStream'));
            
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
            
            % get noisy measurement
            d = model.measure();
            R = model.measureCov();
            
            % sample from measurement error
            if strcmp(this.opts.sampling, 'advanced') && this.opts.samplingFactor > 1.0
                E = model.measureErr(this.opts.samplingFactor*N, this.measurementEnsembleRng);
                E = tools.sampling.condenseSample(E, N);
            else
                E = model.measureErr(N, this.measurementEnsembleRng);
            end
            
            
            if this.opts.samplingOrder >= 1
                % OPTIONAL: first order correct sampling for measurement anomalies/deviations/...
                E = tools.sampling.removeBias(E);
                
                if this.opts.samplingOrder >= 2
                    % OPTIONAL: second order correct sampling for measurement anomalies/deviations/...
                    E = bsxfun(@times, E, sqrt(diag(R)) ./sqrt(var(E,0,2)));
                end
            end
            if strcmp(this.opts.sampleBias, 'remove')
                % OPTIONAL: Remove sampling bias
                E = tools.sampling.removeBias(E);
            end
            
            % create measurement ensemble
            D = bsxfun(@plus,E,d);
            
            % compute ensemble representation of measurement error
            % covariance matrix
%             Re = (1/(N-1)) * (E * E');
            
            % obtain measurements from the ensemble
            h = model.measureOp();
            HA = h(A);
            
            % compute ensemble of innovation vectors
            Dprime = D - HA;
                        
            % substract effect mean from effects on ensemble
            Aprime  =  bsxfun(@minus, A, mean(A, 2));
            HAprime = bsxfun(@minus, HA, mean(HA, 2));
            
            if strcmp(this.opts.measurementCov, 'exact')
                % exact representation of the measurement error covariance
                % matrix
                C = HAprime*HAprime' + (N-1)*R;
            elseif strcmp(this.opts.measurementCov, 'low-rank')
                % low-rank approximation to the measurement error
                % covariance (see Evensen2004, eq. (12) & (14) and
                % Evensen2009a for a discussion)
                C = HAprime*HAprime' + E*E';
            end
            
            if strcmp(this.opts.inflation, 'adaptive')
                % prepare random matrix B for adaptive scaling
                B = randn(floor(this.opts.gamma*M), N);
                B = tools.sampling.removeBias(B);
                for i = 1:floor(this.opts.gamma*M) % scale variance to be exactly 1 in each row
                    B(i,:) = B(i,:) / std(B(i,:));
                end
            end
            
            % Now we have everything ready to compute the real update
            covariance = Aprime * HAprime';
                        
            if strcmp(this.opts.solving, 'MATLAB')
                % Variant 1: Solve using MATLAB backslash operator
                if strcmp(this.opts.inflation, 'adaptive')
                    B = B + Aprime * HAprime' * (C \ Dprime);
                    this.currentRho = mean(std(B,[],2)); % compute inflation factor as the mean over all row standard deviations
                    m = mean(A, 2);
                    A = this.currentRho .* (A - repmat(m, 1, size(A,2))) + repmat(m, 1, size(A,2));
                elseif strcmp(this.opts.inflation, 'fixed') && this.opts.fixedRho > 1.0
                    % directly scale the ensemble perturbation variance
                    m = mean(A, 2);
                    A = this.opts.fixedRho .* (A - repmat(m, 1, size(A,2))) + repmat(m, 1, size(A,2));
                end
                
                A = A + covariance * (C \ Dprime);
                
            elseif strcmp(this.opts.solving, 'eigen')
                % Variant 2: Solve using (limited) eigenvalue decomposition
                % like described in Evensen2003, eqn. (55,56).
                % Important: different types of measurements must be
                % assimilated sequentially, or scaled.
                
                [V, D] = eig(C);
                % [V, D] = eigs(HAprime*HAprime' + E*E',[], min(N,length(d)));
                % eigs is a little bit slower than eig...
                retain_eigs = max(N,length(d));
                weights = V(1:retain_eigs,1:retain_eigs) * diag(1./diag(D(1:retain_eigs,1:retain_eigs))) * V(1:retain_eigs,1:retain_eigs)';
                
                if strcmp(this.opts.inflation, 'adaptive')
                    B = B + Aprime * HAprime' * weights * Dprime;
                    this.currentRho = mean(std(B,[],2)); % compute inflation factor as the mean over all row standard deviations
                    m = mean(A, 2);
                    A = this.currentRho .* (A - repmat(m, 1, size(A,2))) + repmat(m, 1, size(A,2));
                elseif strcmp(this.opts.inflation, 'fixed') && this.opts.fixedRho > 1.0
                    % directly scale the ensemble perturbation variance
                    m = mean(A, 2);
                    A = this.opts.fixedRho .* (A - repmat(m, 1, size(A,2))) + repmat(m, 1, size(A,2));
                end
                
                A = A + covariance * weights * Dprime;
            end
            
            representation.setEnsemble(A, sample);
        end
        
        function stat = getStatistics(this)
            if strcmp(this.opts.inflation, 'adaptive')
                stat.rho = this.currentRho;
            elseif strcmp(this.opts.inflation, 'fixed')
                stat.rho = this.opts.fixedRho;
            else
                stat = struct([]);
            end
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

