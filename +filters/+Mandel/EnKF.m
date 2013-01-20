classdef EnKF < filters.Filter
    %ENKF Standard ensemble Kalman filter (EnKF) according to Mandel2006a
    
    properties
        measurementEnsembleRng;
    end
    
    methods
        function this = EnKF(varargin)
            this.measurementEnsembleRng = RandStream('mlfg6331_64');
        end
        
        function this = update(this, model, representation)
            % transform ensemble to ensemble matrix
            [X sample] = representation.getEnsemble();
            d = model.measure();
            R = model.measureCov();
            
            M = size(X, 1); %#ok<NASGU>
            N = size(X, 2);
            h = model.measureOp();
            HX = h(X);
            
            % compute effect mean
            EHX = mean(HX,2);
            
            % compute ensemble mean
            EX = mean(X, 2);
            
            % substract effect mean from effects on ensemble
            A = bsxfun(@minus, X, EX);
            HA = bsxfun(@minus, HX, EHX);
            
            
            E = model.measureErr(N, this.measurementEnsembleRng);
            E = tools.sampling.removeBias(E);
            E = bsxfun(@times, E, sqrt(diag(R)) ./sqrt(var(E,0,2)));
            
            % create measurement ensemble: add mean to E
            D = bsxfun(@plus, E, d);
            
            % perform the Ensemble Kalman Filter update step
            Y = D - HX;
            P = R + (1/(N-1)) * (HA * HA');
            
            % the following solution should be obtained by a Cholesky
            % decomposition of P
            M = P \ Y;
            representation.setEnsemble(X + ((1 / (N-1)) * A * HA') * M, sample);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

