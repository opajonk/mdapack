classdef PF < filters.Filter
    %PF Particle Filter
    
    properties
        sampling;
        samplingOrder;
        samplingFactor;
        denoise;
        filterOrder;
    end
    
    methods
        function this = PF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('sampling', 'default', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            
            p.addOptional('samplingFactor', 1.0, @(x)x>=1.0 );
            p.addOptional('denoise', false, @(x) isa(x, 'logical'));
            p.addParamValue('samplingOrder', 2, @(x) any(find([0,1,2] == x)));
            
            p.parse(varargin{:});
            
            this.sampling = p.Results.sampling;
            this.samplingFactor = p.Results.samplingFactor;
            this.samplingOrder = p.Results.samplingOrder;
            this.denoise = p.Results.denoise;
        end
        
        function this = update(this, model, representation)
            [coefficients_prior basis] = representation.getPCE();
            [X sample] = representation.getEnsemble();
            N = size(X, 2);
            M = size(X, 1);
            
            % prepare noisy measurement of the truth to assimilate
            d = model.measure();
            
            % simulated measurements of the ensemble
            h = model.measureOp();
            HX = h(X);
            HA = bsxfun(@minus, HX, mean(HX,2));
            R = model.measureCov();
            
            mCov = 1/(N-1)*(HA*HA');
            % compute Gaussian likelihood
            diff = bsxfun(@minus, HX, d);
            weights = zeros(1,N);
            
            E = model.measureErr(N);
            if this.samplingOrder >= 1
                % OPTIONAL: first order correct sampling for measurement anomalies/deviations/...
                E = tools.sampling.removeBias(E);
                
                if this.samplingOrder >= 2
                    % OPTIONAL: second order correct sampling for measurement anomalies/deviations/...
                    E = bsxfun(@times, E, sqrt(diag(R))./sqrt(var(E,0,2)));
                end
            end
            [U S ~] = svd(HA + E,'econ');
            iR=U * diag(1./diag(S*S')*(N-1)) * U';
            
            for i = 1:N
                % FIXME: replace inv by pseudoinverse!
                weights(i) = exp(-0.5 * diff(:,i)'*inv(R)*diff(:,i));
                %                 weights(i) = exp(-0.5 * diff(:,i)'*inv(R + mCov)*diff(:,i));
            end
            
            % normalize
            weights = weights ./ sum(weights);
            ess = this.effectiveSampleSize(weights);
            
            coefficients_post = tools.create_pce(M, representation.pceOrder, X, sample, weights);
            
%             if (ess < 0.5)
%                 % resample
%                 %                 index = this.resampleStratified(weights);
%                 index = this.resampleResidual(weights);
%                 X(:,:)=X(:,index);
%                 sample(:,:)=sample(:,index);
%                 weights=ones(1,N)./N;
%             end
            
%             representation.setEnsemble(X, sample, weights);
            representation.setPCE(coefficients_post, basis);
        end
        
        function ess = effectiveSampleSize(this, w)
            M = size(w,2);
            cv = sum( ( w.*M - 1 ).^2 ) / M;
            ess = 1 / ( 1 + cv );
        end
        
        function index = resampleStratified(this, w)
            N = length(w);
            Q = cumsum(w);
            
            for i=1:N,
                T(i) = rand(1,1)/N + (i-1)/N;
            end
            T(N+1) = 1;
            
            i=1;
            j=1;
            
            while (i<=N),
                if (T(i)<Q(j)),
                    index(i)=j;
                    i=i+1;
                else
                    j=j+1;
                end
            end
        end
        
        function [ indx ] = resampleMultinomial(this, w )
            
            M = length(w);
            Q = cumsum(w);
            Q(M)=1; % Just in case...
            
            i=1;
            while (i<=M),
                sampl = rand(1,1);  % (0,1]
                j=1;
                while (Q(j)<sampl),
                    j=j+1;
                end;
                indx(i)=j;
                i=i+1;
            end
        end
        
        
        function [ indx ] = resampleSystematic(this, w )
            
            N = length(w);
            Q = cumsum(w);
            
            T = linspace(0,1-1/N,N) + rand(1)/N;
            T(N+1) = 1;
            
            i=1;
            j=1;
            
            while (i<=N),
                if (T(i)<Q(j)),
                    indx(i)=j;
                    i=i+1;
                else
                    j=j+1;
                end
            end
        end
        
        
        function [ indx ] = resampleResidual(this, w )
            
            M = length(w);
            
            % "Repetition counts" (plus the random part, later on):
            Ns = floor(M .* w);
            
            % The "remainder" or "residual" count:
            R = sum( Ns );
            
            % The number of particles which will be drawn stocastically:
            M_rdn = M-R;
            
            % The modified weights:
            Ws = (M .* w - floor(M .* w))/M_rdn;
            
            % Draw the deterministic part:
            % ---------------------------------------------------
            i=1;
            for j=1:M,
                for k=1:Ns(j),
                    indx(i)=j;
                    i = i +1;
                end
            end;
            
            % And now draw the stocastic (Multinomial) part:
            % ---------------------------------------------------
            Q = cumsum(Ws);
            Q(M)=1; % Just in case...
            
            while (i<=M),
                sampl = rand(1,1);  % (0,1]
                j=1;
                while (Q(j)<sampl),
                    j=j+1;
                end;
                indx(i)=j;
                i=i+1;
            end
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

