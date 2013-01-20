classdef LEnKF < filters.Filter
    %LENKF Localized ensemble Kalman filter (LEnKF)
    
    properties
        sampling;
        samplingOrder;
        samplingFactor;
        localizationFunction;
        denoise;
        gridPointTransform;
        measurementLocalization;
        measurementRng; % the PRNG used to generate the measurement ensemble
    end
    
    methods
        function this = LEnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('sampling', 'default', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            p.addOptional('samplingFactor', 1.0, @(x)x>=1.0 );
            p.addOptional('samplingOrder', 2, @(x) any(find([0,1,2] == x)));
            p.addOptional('localizationFunction',...
                @(L, dl, N) exp(-1.0 .* (L./dl).^2),...
                @(f) isa(f, 'function_handle'));
            p.addOptional('measurementLocalization', true, @(x) isa(x, 'logical'));
            p.addOptional('denoise', false, @(x) isa(x, 'logical'));
            p.addOptional('gridPointTransform', false, @(x) isa(x, 'logical'));
            p.addOptional('measurementRng', RandStream('mlfg6331_64'), @(x) isa(x,'RandStream'));
            
            p.parse(varargin{:});
            
            this.sampling = p.Results.sampling;
            this.samplingFactor = p.Results.samplingFactor;
            this.samplingOrder = p.Results.samplingOrder;
            this.localizationFunction = p.Results.localizationFunction;
            this.denoise = p.Results.denoise;
            this.gridPointTransform = p.Results.gridPointTransform;
            this.measurementLocalization = p.Results.measurementLocalization;
            this.measurementRng = p.Results.measurementRng;
        end
        
        function this = update(this, model, representation)
            % transform ensemble to ensemble matrix
            [X sample] = representation.getEnsemble();
            % check validity of ensemble
            assert(isempty(find(~logical(isfinite(X)), 1)), 'the INPUT ensemble is invalid as it contains NaN or Inf values - we have to abort');
            
            % Obtain some helper variables, the measurements,
            % and the grid space localizer from the model.
            h = model.measureOp();
            N = size(X, 2);
            L = this.localizationFunction(model.distanceMatrix(), model.decorrelationLength(), N);
            HX = h(X);
            
            % This method locates the maximum of the measurement operator and
            % "assigns" the measurement location to that maximum in
            % the model domain. In case of point measurements this is the
            % same as directly using the measurement operator.
            % The resulting H we use to construct a localizer for
            % the measurements, mL. Another possibility could be to use the
            % "support". See also Fertig2007a, appendix A&B.
            mpos = tools.computeMPos(model.measurementDimension, model.deterministicDimension, model.measureOp());
            H = zeros(model.measurementDimension, model.deterministicDimension);
            for i = 1:model.measurementDimension
                H(i,mpos(i)) = 1;
            end
            mL = H*L;
            
            if this.gridPointTransform == true
                % Build "gridpoint transform" which represents the smallest grid
                % point heterogeneities. Such we obtain a correlation matrix
                % later on, e.g. when computing A*HA', instead of a covariance
                % matrix. It also "removes the dimension" from grid-point
                % and measurement space, obliviating the need for
                % "normalisation" of measurements for proper SVD-based
                % inversion.
                % See also Deckmyn2005, p. 1283.
                A = bsxfun(@minus, X, mean(X, 2));
                HA = bsxfun(@minus, HX, mean(HX, 2));
                B = spdiags(std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
                iB = spdiags(1./std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
                iHB = spdiags(1./std(HA,0,2),0,model.measurementDimension,model.measurementDimension);
                
                X = iB*X;
                HX = iHB*HX;
            end
            
            % optional denoising step of the ensemble using Wavelets
            if this.denoise == true
                for i = 1:size(X,2)
                    X(:,i) = wden(X(:,i),'heursure','s','one',3,'sym8');
                end
            end

            % obtain the noisy measurement
            if this.gridPointTransform == true
                d = iHB*model.measure();
            else
                d = model.measure();
            end
            
            % obtain the measurement covariance
            if this.gridPointTransform == true
                R = iHB*model.measureCov()*iHB';
            else
                R = model.measureCov();
            end
            
            % compute ensemble deviations/anomalies/whatever-they-call-this...
            A = bsxfun(@minus, X, mean(X, 2));
            HA = bsxfun(@minus, HX, mean(HX,2));
            
            % sample from measurement error
            if strcmp(this.sampling, 'advanced') && this.samplingFactor > 1.0
                % Advanced sampling
                if this.gridPointTransform == true
                    E = iHB*model.measureErr(this.samplingFactor*N, this.measurementRng);
                else
                    E = model.measureErr(this.samplingFactor*N, this.measurementRng);
                end
                E = tools.sampling.condenseSample(E, N);
            else
                % Default sampling
                if this.gridPointTransform == true
                    E = iHB*model.measureErr(N,this.measurementRng);
                else
                    E = model.measureErr(N,this.measurementRng);
                end
            end
            
            if this.samplingOrder >= 1
                % OPTIONAL: first order correct sampling for measurement anomalies/deviations/...
                E = tools.sampling.removeBias(E);
                
                if this.samplingOrder >= 2
                    % OPTIONAL: second order correct sampling for measurement anomalies/deviations/...
                    E = bsxfun(@times, E, sqrt(diag(R))./sqrt(var(E,0,2)));
                end
            end
            
            % create measurement ensemble
            D = bsxfun(@plus, E, d);
            
            % Perform the Ensemble Kalman Filter update step using
            % pseudo-inversion. This is necessary to make the covariance localization
            % possible.
            if this.measurementLocalization == true
                C = (HA * HA');
                C = (mL) .* C;
                P = (N-1).*R + C;
            else
                P = (N-1).*R + (((HA * HA')));
            end
            
            C = A * HA';
            C = L.* C;
            K = C * pinv(P);
            
            X = X + K * (D - HX);
            
            % transform back into grid point space if necessary
            if this.gridPointTransform == true
                X = B*X;
            end
            
            % check validity of the output ensemble
            assert(isempty(find(~logical(isfinite(X)), 1)), 'the OUTPUT ensemble is invalid as it contains NaN or Inf values - we have to abort');
            
            representation.setEnsemble(X, sample);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

