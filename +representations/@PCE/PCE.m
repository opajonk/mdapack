classdef PCE < representations.Representation
    %PCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coefficients;
    end
    
    properties (SetAccess = private)
        basis;
        mTensor;
        pceOrder;
        sampleSize;
        sample;
        samplingOrder;
        multTensor;
        multFactorial;
        rndStream;
    end
    
    methods
        function this = PCE(model, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addRequired('model', @(x) isa(x, 'models.Model'));
            p.addParamValue('pceOrder', 1, @(x)x>=0)
            p.addParamValue('sampleSize', 1000, @(x)x>=10)
            p.addParamValue('samplingOrder', 2, @(x) any(find([0,1,2] == x)));
            p.addParamValue('rndStream', RandStream('mt19937ar','Seed',tools.shuffleSeed), @(x) isa(x,'RandStream'));
            
            p.parse(model, varargin{:});
            
            this.pceOrder = p.Results.pceOrder;
            this.sampleSize = p.Results.sampleSize;
            this.samplingOrder = p.Results.samplingOrder;
            this.rndStream = p.Results.rndStream;
            
            this.basis = multiindex(model.stochasticDimension, this.pceOrder);
            this.multTensor = hermite_triple_fast(this.basis,this.basis,this.basis,'algorithm','default');
            this.multFactorial = multiindex_factorial(this.basis)';
            model.createInitialPCE(this);
            
            this.sample = randn(this.rndStream, size(this.basis, 2), this.sampleSize);
            
            if this.samplingOrder >= 1
                this.sample = tools.sampling.removeBias(this.sample);
                
                if this.samplingOrder >= 2
                    this.sample = bsxfun(@times, this.sample, 1.0./sqrt(var(this.sample,0,2)));
                end
            end
        end
        
        function [coefficients basis multFactorial multTensor] = getPCE(this)
            coefficients = this.coefficients;
            basis = this.basis;
            multTensor = this.multTensor;
            multFactorial = this.multFactorial;
        end
        
        function [this] = setPCE(this, coefficients, basis)
            this.coefficients = coefficients;
            this.basis = basis;
            this.multTensor = hermite_triple_fast(this.basis,this.basis,this.basis);
            this.multFactorial = multiindex_factorial(this.basis)';
        end
        
        function [ensemble sample] = getEnsemble(this)
            ensemble = pce_evaluate(this.coefficients, this.basis, this.sample);
            sample = this.sample;
        end
        
        function [this] = setEnsemble(this, ensemble, sample, weights)
            N = size(ensemble, 2);
            
            % default weights: all equal
            if nargin < 4
                weights = ones(1, N) ./ N;
            end
            this.coefficients = tools.create_pce(size(ensemble, 1), this.pceOrder, ensemble, sample, weights);
        end
        
        function ret = mean(this)
            ret = this.coefficients(:,1);
        end
        
        function ret = summary(this)
%             sample = randn(this.rndStream, size(this.basis, 2), this.sampleSize);
%             
%             if this.samplingOrder >= 1
%                 sample = tools.sampling.removeBias(sample);
%                 
%                 if this.samplingOrder >= 2
%                     sample = bsxfun(@times, sample, 1.0./sqrt(var(sample,0,2)));
%                 end
%             end
            ensemble = pce_evaluate(this.coefficients, this.basis, this.sample);
            ret(:,1) = tools.percentile(ensemble', .025);
            ret(:,2) = tools.percentile(ensemble', .25 );
            ret(:,3) = tools.percentile(ensemble', .50 );
            ret(:,4) = tools.percentile(ensemble', .75 );
            ret(:,5) = tools.percentile(ensemble', .975);
        end
        
        function ret = var(this)
            % direct computation, not by integration or something
            ret=this.coefficients(:,2:end).^2*this.multFactorial(:,2:end)';
        end
        
        function ret = std(this)
            ret = sqrt(this.var());
        end
        
        function ret = rmse(this, truth)
            factor = this.coefficients;
            factor(:,1) = factor(:,1) - truth;
            
            ret=sqrt(mean(factor.^2*this.multFactorial(:,1:end)'));
        end
        
        function ret = skew(this)
            [~,~,ret,~] = pce_moments(this.coefficients, this.basis);
        end
        
        function ret = kurt(this)
            [~,~,~,ret] = pce_moments(this.coefficients, this.basis);
        end
        
        function ret = moment(this, order)
            assert(order > 0 && order < 5);
            [MEAN,VAR,SKEW,KURT] = pce_moments(this.coefficients, this.basis);
            
            switch(order)
                case 1
                    ret = MEAN;
                case 2
                    ret = VAR;
                case 3
                    ret = SKEW;
                case 4
                    ret = KURT;
            end
        end
    end
end

