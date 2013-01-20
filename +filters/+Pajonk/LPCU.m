classdef LPCU < filters.Filter
    %LPCU Linear Polynomial Chaos Update
    
    properties
        opts;
    end
    
    methods
        function this = LPCU(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('square_root_scheme_enabled', true, @(x) isa(x, 'logical'));
            p.addOptional('square_root_scheme_symmetric', true, @(x) isa(x, 'logical'));
            p.addOptional('Zeng2010', false, @(x) isa(x, 'logical'));
            p.parse(varargin{:});
            
            this.opts = p.Results;
        end
        
        function this = update(this, model, representation)
            [X I_X mF_X] = representation.getPCE();
            
            h = model.measureOp();
            d = model.measure();
            R = model.measureCov();
            
            N = size(X, 2);
            d_dim = length(d);
            
            HX = h(representation);
            
            % Multiply by the square root of the gramian so that the
            % variance scales linearly with the PCE coefficients. Then
            % nHX and nX can be used for cheap computation of covariance
            % matrices, as well as in the square root scheme.
            nHX = bsxfun(@times, HX(:,2:end), sqrt(mF_X(2:end)));
            nX = bsxfun(@times, X(:,2:end), sqrt(mF_X(2:end)));
            
            P = nHX * nHX' + R;
            C_XHX = nX * nHX';
            K = C_XHX * pinv(P);
            
            if (this.opts.square_root_scheme_enabled)
                Xup = zeros(size(X));
                
                % First update the mean as usual
                Xup(:,1) = X(:,1) + K * (d - HX(:,1)); 
                
                % Then update the variance by updating the PCE coefficients
                % with a square root scheme.
                [Z, LAMBDA] = eig(P);
                nHX2 = diag(1./diag(sqrt(LAMBDA))) * Z' * nHX;
                [~, SIGMA, V] = svd(nHX2);
                
                % There are many different implementation possibilities of
                % the square root scheme.
                if this.opts.square_root_scheme_symmetric
                    % This is the "symmetric" square root scheme, which
                    % re-distributes the variance across the full PCE.
                    nX2 = nX * V * diag(sqrt(diag(eye(N-1) - SIGMA' * SIGMA))) * V';
                else
                    % This is probably the simplest one. But it always
                    % gives a multivariate Gaussian posterior (at least in
                    % the case when dim(obs)==dim(state))
                    nX2 = nX *  V * diag(sqrt(diag(eye(N-1) - SIGMA' * SIGMA)));
                end
                
                % Fix: the measurement RVs may be fully orthogonal to one
                % or more prior RVs, causing the posterior to have 0 variance
                % due to the limited update. However, prior = posterior is
                % correct in that case! So set this:
                zeroVarRVs = find(sum(nX2.^2,2) == 0);
                if (~isempty(zeroVarRVs))
                    nX2(zeroVarRVs,:) = nX(zeroVarRVs,:);
                end
                
                % Divide out the square root of the gramian, so that the
                % result is a valid PCE
                Xup(:,2:end) = bsxfun(@rdivide, nX2, sqrt(mF_X(2:end)));
            else
                D = zeros(size(HX));
                D(:,1) = d;
                
                % Zeng2010 use Dirac delta for the PCE terms of the evidence in the
                % difference D-HX.
                if (norm(full(R)) > 0 &&  ~this.opts.Zeng2010)
                    D(:,2:d_dim+1) = chol(R);
                end

                innovations = K * (D-HX);
                
                % Is that actually necessary?
                innovations = bsxfun(@rdivide, innovations, sqrt(mF_X));
                
                if size(innovations) ~= size(X)
                    innovations = innovations';
                end

                Xup = X + innovations;
            end
            
            representation.setPCE(Xup, I_X);
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

