classdef EnSRF < filters.Filter
    %EnSRF Square root implementation of the EnKF (e.g. Evensen2004,
    %Evensen2009a)
    
    properties
        opts;
    end
    
    methods
        function this = EnSRF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('squareRoot', 'symmetric', ...
                @(x)any(strcmpi(x,{'symmetric','standard'})));
            p.addOptional('randomRotation', false, ...
                @(x)islogical(x));
            
            p.parse(varargin{:});
            
            this.opts = p.Results;
        end
        
        function this = update(this, model, representation)
            % preparations
            [A sample] = representation.getEnsemble();
            M = size(A, 1);
            N = size(A, 2);
            
            % get noisy measurement
            d = model.measure();
            
            % compute ensemble representation of measurement error
            % covariance matrix
            % Re = (1/(N-1)) * (E * E');
            R = model.measureCov();
            
            % obtain measurements from the ensemble
            h = model.measureOp();
            HA = h(representation);
            
            % compute measurement mean
            meanHA = mean(HA,2);
            
            % compute ensemble mean
            meanA = mean(A, 2);
            
            % substract effect mean from effects on ensemble
            Aprime  =  bsxfun(@minus, A, mean(A, 2));
            HAprime = bsxfun(@minus, HA, mean(HA, 2));
            
            % now we have ready: A, HA, A', HA', d, R
            C = HAprime*HAprime' + (N-1)*R;
            
            [Z, LAMBDA] = eig(C);
            X2 = diag(1./diag(sqrt(LAMBDA))) * Z' * HAprime;
            [~, SIGMA, V] = svd(X2);
            
            
            if strcmp(this.opts.squareRoot, 'symmetric')
                % this is the "symmetric" square root implementation (e.g.
                % Evensen2009a, eq. (75) and Sakov2008, Livings2008
                Aprime_analyzed = Aprime * V * diag(sqrt(diag(eye(N) - SIGMA' * SIGMA))) * V';
            else
                % this is the "standard" or "Evensen2004" implementation
                Aprime_analyzed = Aprime * V * diag(sqrt(diag(eye(N) - SIGMA' * SIGMA)));
            end
            
            if (this.opts.randomRotation)
                % Apply a random rotation. This can be seen as a
                % "resampling from posterior" procedure. See Evensen2009a,
                % eq. (76) and Evensen2004, eq, (34) and Evensen2009, p.
                % 200 and eq. (13.12)
                [~,~,THETA] = svd(randn(N),'econ');
                Aprime_analyzed = Aprime_analyzed * THETA';
            end
            
            meanA_analyzed = meanA + Aprime * HAprime' * (Z * diag(1./diag(LAMBDA)) * Z')*(d - meanHA);
            
            A = bsxfun(@plus, Aprime_analyzed, meanA_analyzed);
            
            representation.setEnsemble(A, sample);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

