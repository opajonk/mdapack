classdef MHF < filters.Filter
    %MHF Metropolis-Hastings-based filter
    
    properties
        opts;
        
        proposalRng;
    end
    
    methods
        function this = MHF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('nsamples', 10000, @(x)x>=100);
            p.addOptional('proposalStdDev', 0.1, @(x)x>0.0);
            p.addOptional('proposalRng', RandStream('mlfg6331_64'), @(x) isa(x,'RandStream'));
            
            p.parse(varargin{:});
            
            this.opts = p.Results;
            
            % extract the PRNG from the options - it is not an option, so
            % remove it from there
            this.proposalRng = this.opts.proposalRng;
            this.opts = rmfield(this.opts, 'proposalRng');
        end
        
        function this = update(this, model, representation)
            [coefficients_prior basis] = representation.getPCE();
            M = model.stochasticDimension;
                        
            % Prepare data which needs to be acquired from the model.
            d = model.measure();
            R = model.measureCov();
            h = model.measureOp();
            
            % Prepare the necessary PDFs for Metropolis-Hastings.
            % Here we use a Gaussian centered at the curret position as
            % proposal RNG and distribution.
            logtargetpdf = @(x) (this.posteriorPdf(x, coefficients_prior, basis, d, h, inv(R)));
            logproposalpdf = @(x,y) -0.5 *((x-y)*inv(this.opts.proposalStdDev.^2)*(x-y)'); %#ok<MINV>
            logproposalrnd = @(x) x + this.opts.proposalStdDev*this.proposalRng.randn(size(x));
            
            % Use fminsearch to find a good start value - otherwise it may
            % happen that MHSAMPLE bails out with an error.
            start = zeros(1,model.stochasticDimension);

            % Construct Metropolis-Hastings-chain - these are samples from the
            % Bayesian posterior pdf
            [x accept] = mhsample(start,...
                this.opts.nsamples,...
                'logpdf', logtargetpdf,...
                'logproppdf', logproposalpdf,...
                'proprnd', logproposalrnd,...
                'symmetric', 1,...
                'burnin', floor(this.opts.nsamples*0.2), ... % 20% burnin
                'thin', 3); % thin chain to lower repeated samples
            
            % Just check if the acceptance rate is OK.
            % "The asymptotically optimal acceptance rate is 0.234 under
            % quite general conditions." (see [Gelman1997])
            
            if (accept > 0.95)
                fprintf(1,'INFO: Very high acceptance rate of MCMC (%f) - recommended is  ~0.234\n', accept);
            elseif (accept < 0.05)
                fprintf(1,'INFO: Very low acceptance rate of MCMC (%f) - recommended is  ~0.234\n', accept);
%             else
%                 fprintf(1,'INFO: Acceptance rate of MCMC is %f\n', accept);
            end
            
            % Compute samples from the posterior in "model space" by
            % evaluating the PCE.
            ynew = pce_evaluate(coefficients_prior, basis, x');
            
            % The samples x of the posterior are samples of the thetas (the
            % standard normal RVs which go into the PCE) - which are of
            % course no longer standard normal in the posterior!
            % Construct the Nataf transformation to map them to
            % standard normal, uncorrelated (independent) RVs,
            % then construct new PCE from those transformed samples xnew
            % and above "samples in model space" ynew - that is the
            % posterior PCE which we will feed back to the model.
            
            N = size(x,1);
            
            % Perform Nataf transform (see e.g. [Lebrun2009a])
            xnew = tools.nataf_transform(x);
%             xnew = x';

%             Perform empirical Rosenblatt
            for k = 1:M
                [hcdf] = tools.buildhcdf(xnew(1:k,:));
                for i=1:N
                    xnew(k,i) = tools.evalhcdf(hcdf, xnew(1:k,i));
                end
            end
            
            xnew = tools.nataf_transform(xnew');
            
            % Perform empirical Rosenblatt in reverse direction
            for k = M:1
                [hcdf] = tools.buildhcdf(xnew(M:k,:));
                for i=1:N
                    xnew(k,i) = tools.evalhcdf(hcdf, xnew(M:k,i));
                end
            end
            
            % Perform Nataf transform 
            xnew = tools.nataf_transform(xnew');
            
            % Perform inverse cdf transform
            for i=1:M
                [~,~,xmesh,cdf] = tools.kde(xnew(i,:));
                
                % HACK(?): Make sure the values of the density estimate are
                % valid - otherwise xnew may contain NaNs
                % Is there a better solution for this?
                cdf(cdf<=0.0) = eps;
                cdf(cdf>=1.0) = 1-eps;
                
                % interpolate using piecewise cubic Hermite polynomials
                % these preserve monotonicity (in the way they are setup in MATLAB)
                xnew(i,:) = pchip(xmesh,cdf,xnew(i,:));
            end
            
            % Perform transform to Gaussian marginals
            for k = 1:M
                xnew(k,:) = norminv(xnew(k,:),0,1);
            end
            
            % Construct posterior PCE from these samples
            coefficients_post = tools.create_pce(M, representation.pceOrder, ynew, xnew);
%             coefficients_post(4,:) = coefficients_prior(4,:);
            
            % TEST
%             samples = randn(size(xnew));
%             ytest = pce_evaluate(coefficients_prior, basis, samples);
            
            representation.setPCE(coefficients_post, basis);
        end
        
        
        function [l] = posteriorPdf(this, x, coefficients, basis, data, h, iR)
            y = pce_evaluate(coefficients, basis, x');
            diff = data - h(y);
            %             l = exp(-0.5 * (diff'*iR*diff)); % likelihood...
            %             l = l * exp(-0.5 *(x*x')); % * prior: N(0,1) == prior PCE
            l = -0.5 * (diff'*iR*diff); % likelihood...
            l = l + -0.5 *(x*x'); % * prior: N(0,1) == prior PCE
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

