classdef JDEnKF < filters.Filter
    %JDEnKF Joint diagonalisation ensemble Kalman filter (JDEnKF)
    
    properties
        sampling;
        samplingOrder;
        samplingFactor;
        localizationFunction;
        denoise;
        measurementRng; % the PRNG used to generate the measurement ensemble
    end
    
    methods
        function this = JDEnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('sampling', 'default', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            
            p.addOptional('samplingFactor', 1.0, @(x)x>=1.0 );
            p.addOptional('denoise', false, @(x) isa(x, 'logical'));            
            p.addOptional('samplingOrder', 2, @(x) any(find([0,1,2] == x)));
            p.addOptional('measurementRng', RandStream('mlfg6331_64'), @(x) isa(x,'RandStream'));
            p.addOptional('localizationFunction',...
                @(L, dl, N) exp(-1.0 .* (L./dl).^2),...
                @(f) isa(f, 'function_handle'));
            
            p.parse(varargin{:});
            
            this.sampling = p.Results.sampling;
            this.samplingFactor = p.Results.samplingFactor;
            this.samplingOrder = p.Results.samplingOrder;
            this.localizationFunction = p.Results.localizationFunction;
            this.denoise = p.Results.denoise;
            this.measurementRng = p.Results.measurementRng;
        end
        
        function this = update(this, model, representation)
            % transform ensemble to ensemble matrix
            [X sample] = representation.getEnsemble();
            
            % optional denoising step of the ensemble using Wavelets
            if this.denoise == true
                for i = 1:size(X,2)
                    X(:,i) = wden(X(:,i),'heursure','s','one',3,'sym8');
                end
            end
            
            % prepare noisy measurement of the truth to assimilate
            d = model.measure();
            h = model.measureOp();
            % simulated measurements of the ensemble
            HX = h(representation);
            H = h(eye(model.deterministicDimension));
            
            R = model.measureCov();
            N = size(X, 2);
            
            % compute ensemble deviations/anomalies/whatever-they-call-this...
            A = bsxfun(@minus, X, mean(X, 2));
            HA = bsxfun(@minus, HX, mean(HX,2));
                        
            % this operator creates is a transformation to a correlation
            % matrix with diagonal normalized to 1.0
            % (unused now)
%             B = spdiags(std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
%             iB = spdiags(1./std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
%             HB = spdiags(std(HA,0,2),0,model.measurementDimension,model.measurementDimension);
%             iHB = spdiags(1./std(HA,0,2),0,model.measurementDimension,model.measurementDimension);

            
            % sample from measurement error
            if strcmp(this.sampling, 'advanced') && this.samplingFactor > 1.0
                % Advanced sampling
                E = model.measureErr(this.samplingFactor*N, this.measurementRng);
                E = tools.sampling.condenseSample(E, N);
            else
                % Default sampling
                E = model.measureErr(N, this.measurementRng);
            end
            
            if this.samplingOrder >= 1
                % OPTIONAL: first order correct sampling for measurement anomalies/deviations/...
                E = tools.sampling.removeBias(E);
                
                if this.samplingOrder >= 2
                    % OPTIONAL: second order correct sampling for measurement anomalies/deviations/...
                    E = bsxfun(@times, E, sqrt(diag(R))./sqrt(var(E,0,2)));
                end
            end
        
            % compute a simultaneous diagonalisation of the measurement
            % operator H*H' and the measurement covariance HA*HA', E*E'
            [W] = tools.joint_diag_r([(R) (HA*HA') (H*H')]);
            
            HAd = diag(W'*(HA*HA')*W);
            Hd = diag(W'*(H*H')*W);
            Ed = diag(W'*(R)*W);
            
            % create measurement ensemble
            D = bsxfun(@plus, E, d);
            
            % join diagonalisation approximate analysis
            X = X + A*HA' * (W*diag(1./(HAd + Ed))*W') * (D - HX);
            
            representation.setEnsemble(X, sample);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

