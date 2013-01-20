classdef KF < filters.Filter
    %KF Original Kalman filter, based on PCE of order 1
    
    properties
        opts;
    end
    
    methods
        function this = KF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
                        
            p.parse(varargin{:});
            
            this.opts = p.Results;
        end
        
        function this = update(this, model, representation)
            [X basis multFactorial multTensor] = representation.getPCE();
            n = size(X, 1);
            h = model.measureOp();
            HX = h(representation);
            
            d = model.measure();
            R = model.measureCov();
            
            %% Kalman filter update
            % compute covariance from PCE
            covHX=double(lbucov(HX,HX,multTensor,multFactorial));
            covXHX=double(lbucov(X,HX,multTensor,multFactorial));
            covXX=double(lbucov(X,X,multTensor,multFactorial));
            
            K = covXHX * pinv(double(covHX)+R);
            
            newX = zeros(size(X));
            
            newX(:,1) = X(:,1) + K * (d - HX(:,1)); % update the mean
            covXX  = (eye(n) - h(K')') * covXX; % update variance
            
            % compute non-mean PCE coefficients of a multivariate Gaussian from the
            % covariance
            pcecoeff = sqrtm(covXX);
            
            newX(:,2:n+1) = pcecoeff;
            
%             test =double(lbucov(newX,newX,multTensor,multFactorial));
%             if (norm(test - covXX) > eps)
%                 disp('this should never happen');
%             end
            
            representation.setPCE(newX, basis);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

