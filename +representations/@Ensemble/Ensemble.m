classdef Ensemble < representations.Representation
    %ENSEMBLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ensembleM; % the current ensemble of state/parameter estimates
        sampleM; % the sample that was used to create the current ensemble (necessary for some algorithms)
        sampling;
        distribution;
        beta;
        sampleSize;
        samplingOrder;
        pceOrder;
        rndStream;
    end
    
    methods (Static = true)
        
    end
    
    methods
        function this = Ensemble(model, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addRequired('model', @(x) isa(x, 'models.Model'));
            
            p.addParamValue('sampleSize', 100, @(x)x>1);
            p.addParamValue('sampling', 'default', @(x)any(strcmpi(x,{'default','advanced'})));
            p.addParamValue('distribution', 'normal', @(x)any(strcmpi(x,{'normal','uniform'})));
            p.addParamValue('beta', 1, @(x)x>=1);
            p.addParamValue('initialEnsembleFileName', '', @(x) isa(x,'char'));
            p.addParamValue('samplingOrder', 2, @(x) any(find([0,1,2] == x)));
            p.addParamValue('pceOrder', 1, @(x)x>=1)
            p.addParamValue('rndStream', RandStream('mt19937ar','Seed', now()), @(x) isa(x,'RandStream'));
            
            p.parse(model, varargin{:});
            
            this.sampleSize = p.Results.sampleSize;
            this.sampling = p.Results.sampling;
            this.beta = p.Results.beta;
            this.distribution = p.Results.distribution;
            this.samplingOrder = p.Results.samplingOrder;
            this.pceOrder = p.Results.pceOrder;
            this.rndStream = p.Results.rndStream;
            
            if ~strcmp(p.Results.initialEnsembleFileName,'')
                % we load a pre-created ensemble matrix from file and extract "sampleSize" members
                in = load(p.Results.initialEnsembleFileName);
                if exist('in.ensemble','var')
                    assert(size(in.ensemble, 2) >= this.sampleSize, 'ensemble is too small'); % matrix must be big enough
                    this.ensembleM = in.ensemble(:,1:this.sampleSize);
                    this.generatingSampleM = in.generatingSample(:,1:this.sampleSize);
                    clear in;
                else
                    error('the supplied file should contain a variable called ensemble');
                end
            else
                if strcmp(this.distribution,'normal')
                    this.sampleM = this.rndStream.randn(model.stochasticDimension, this.sampleSize);
                    if this.samplingOrder >= 1
                        % OPTIONAL: first order correct sampling for measurement anomalies/deviations/...
                        this.sampleM = tools.sampling.removeBias(this.sampleM);
                        
                        if this.samplingOrder >= 2
                            % OPTIONAL: second order correct sampling for measurement anomalies/deviations/...
                            this.sampleM = bsxfun(@times, this.sampleM, 1.0./sqrt(var(this.sampleM,0,2)));
                        end
                    end
                elseif strcmp(this.distribution,'uniform')
                    this.sampleM = rand(this.rndStream, model.stochasticDimension, this.sampleSize);
                end
                
                model.createInitialEnsemble(this);
            end
        end
        
        function ens_size = ens_size(this)
            ens_size = size(this.ensembleM,2);
        end
        
        function ret = mean(this)
            ret = mean(this.ensembleM, 2);
        end
        
        function ret = std(this)
            ret = std(this.ensembleM,[],2);
        end
        
        function ret = var(this)
            ret = var(this.ensembleM,0,2);
        end
        
        function ret = summary(this)
            ret(:,1) = tools.percentile(this.ensembleM', .025);
            ret(:,2) = tools.percentile(this.ensembleM', .25 );
            ret(:,3) = tools.percentile(this.ensembleM', .50 );
            ret(:,4) = tools.percentile(this.ensembleM', .75 );
            ret(:,5) = tools.percentile(this.ensembleM', .975);
        end
        
        function [ensemble sample] = getEnsemble(this)
            ensemble = this.ensembleM;
            sample = this.sampleM;
        end
        
        function [this] = setEnsemble(this, ensemble, sample)
            this.ensembleM = ensemble;
            this.sampleM = sample;
        end
        
        function [coefficients basis multFactorial multTensor] = getPCE(this)
            basis = multiindex(size(this.sampleM, 1), this.pceOrder);
            coefficients = tools.create_pce(size(this.ensembleM, 1), this.pceOrder, this.ensembleM, this.sampleM);
            multTensor = hermite_triple_fast(basis,basis,basis);
            multFactorial = multiindex_factorial(basis)';
        end
        
        function [this] = setPCE(this, coefficients, basis)
            this.ensembleM = pce_evaluate(coefficients, basis, this.sampleM);
        end
        
        function ret = rmse(this, truth)
            factor = this.ensembleM;
            factor = bsxfun(@minus, factor, truth);
            ms = mean(factor.^2,2);
            ret=sqrt(mean(ms));
        end
        
        function ret = skew(this)
            dim = find(size(this.ensembleM) ~= 1, 1, 'last');
            N = size(this.ensembleM, dim);
            ret = this.moment(3);
            ret = ret ./ sqrt(this.var().^3);
            ret = ret .* sqrt((N-1)./N) .* N./(N-2);
        end
        
        function ret = kurt(this)
            dim = find(size(this.ensembleM) ~= 1, 1, 'last');
            N = size(this.ensembleM, dim);
            ret = this.moment(4) ./ this.var().^2;
            ret = ((N+1).* ret - 3.*(N-1)) .* (N-1)./((N-2).*(N-3)) + 3;
            ret = ret - 3; % we follow the convention to substract 3
        end
        
        function ret = moment(this, order)
            assert(order > 0);
            ens = this.ensembleM;
            
            dim = find(size(ens) ~= 1, 1, 'last');
            if isempty(dim), dim = 1; end
            
            meanx = mean(ens,dim);
            
            if (order == 1)
                ret = meanx;
                return;
            end
            
            tile = ones(1,max(ndims(ens),dim));
            tile(dim) = size(ens,dim);
            
            ret = mean((ens - repmat(meanx, tile)).^order,dim);
        end
    end
end

