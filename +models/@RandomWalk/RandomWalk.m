classdef RandomWalk < models.Model
    
    properties (SetAccess = private)
        deterministicDimension;
        stochasticDimension;
        measurementDimension;
        measurementFreq;
        tCount;
        tStep;
        measurementStdDev;
        modelStdDev;
        
        modelNoiseRng;
        evidenceNoiseRng;
        
        h;
    end
    
    properties (Access = private)
        time = 0;
        R;
        mPosData;
        ensMeanData;
        ensStdData;
        summaryData;
        truthData;
        timeScale;
    end
    
    methods (Access = private)        
        function R = computeR(this)
            R = speye(this.measurementDimension) * this.measurementStdDev.^2;
        end
    end
    
    methods (Access = protected)
        function state = createTruth(this)
            state = 123.76721401; % just some random but fixed number
        end
    end
    
    methods
        function this = RandomWalk(varargin)
            p = inputParser;
            
            p.addParamValue('deterministicDimension', 1, @(x)x>1);
            p.addParamValue('measurementDimension', 1, @(x)x>=1);
            p.addParamValue('measurementFreq', 5, @(x)x>0);
            p.addParamValue('tCount', 200, @(x)x>0);
            p.addParamValue('tStep', 1, @(x)x>0);
            p.addParamValue('measurementStdDev', 0.1, @(x)x>0);
            p.addParamValue('modelStdDev', 0.15, @(x)x>0);
            p.addParamValue('measurementOperator', @(x) tools.customMeasurementOperator(x, 1), @(f) isa(f,'function_handle'));
            
            p.addParamValue('modelNoiseRng', RandStream('mt19937ar','seed', 234), @(x) isa(x,'RandStream'));
            p.addParamValue('evidenceNoiseRng', RandStream('mt19937ar','Seed',tools.shuffleSeed), @(x) isa(x,'RandStream'));
            
            p.parse(varargin{:});
            
            this.deterministicDimension = p.Results.deterministicDimension;
            this.measurementDimension = p.Results.measurementDimension;
            this.measurementFreq = p.Results.measurementFreq;
            this.tCount = p.Results.tCount;
            this.tStep = p.Results.tStep;
            this.measurementStdDev = p.Results.measurementStdDev;
            this.modelStdDev = p.Results.modelStdDev;
            this.modelNoiseRng = p.Results.modelNoiseRng;
            this.evidenceNoiseRng = p.Results.evidenceNoiseRng;
            this.h = p.Results.measurementOperator;
            
            assert(this.measurementFreq >= this.tStep, 'The measurement frequency cannot be smaller than the time step!');
            assert(mod(this.measurementFreq, this.tStep) < eps, 'The measurement frequency must be divisible by the time step!');
            
            this.stochasticDimension = this.deterministicDimension;
            this.R = this.computeR();
            
            % create truth case
            this.truth = this.createTruth();
        end
        
        function [coefficients] = createInitialPCE(this, pce) %#ok<STOUT>
            assert(isa(pce, 'representations.PCE'), 'createInitialPCE(): you must pass a PCE to this function');
            assert(size(pce.basis, 2) == this.stochasticDimension, 'createInitialPCE(): wrong stochastic basis')

            pce.coefficients = zeros(size(pce.basis'));
            pce.coefficients(1) = 123.0; % mean
            pce.coefficients(2) = 1.0; % initial std dev
        end
        
        function [] = createInitialEnsemble(this, ens, varargin)
%         function [state] = createInstances(this, sourceSamples)
            assert(isa(ens, 'representations.Ensemble'), 'createInitialEnsemble(): you must pass an Ensemble to this function');
            assert(size(ens.sampleM, 1) == this.stochasticDimension, 'createInitialEnsemble(): wrong stochastic dimension')

            N = size(ens.sampleM, 2);
            
            firstguess = 123.0;
            
            state = firstguess + ens.sampleM;
            
            % rescale initial ensemble for exact variance of
            % 1.0
            m = mean(state, 2);
            
            state = bsxfun(@minus, state, m);
            state = bsxfun(@rdivide, state, std(state, 0, 2));
            state = bsxfun(@plus, state, m);
            ens.ensembleM = state;
        end
        
        function this = step(this, representation, t)
            if nargin<2 || isempty(t)
                t = this.tStep;
            end
            assert(t >= 1);
            
            % advance the truth using the internal random number generator
            % so that the course is deterministic
            this.truth = this.truth + this.modelNoiseRng.randn(1) * this.modelStdDev;
            
            % advance the estimator representation - different implementations may be possible here
            if isa(representation, 'representations.Ensemble')
                sample = randn(size(representation.ensembleM));
                sample = tools.sampling.removeBias(sample);
                sample = bsxfun(@times, sample, this.modelStdDev./sqrt(var(sample,0,2)));
                representation.ensembleM = representation.ensembleM + sample;
            elseif isa(representation, 'representations.PCE')
                % add an independent sample and project back - this amouts
                % to adding stddevs of Gaussians, which goes like sqrt(stddev_a^2 + stddev_b^2)                
                representation.coefficients = [representation.coefficients(1) sqrt(representation.coefficients(2)^2 + this.modelStdDev^2) zeros(length(representation.coefficients) - 2)];
            end            
            
            this.time = this.time + 1;
        end
        
        function [err] = stepErr(this, t)
            assert(false, 'not implemented');
        end
        
        function [hasM] = hasMeasurement(this)
%             if (this.time == 0.0) % this piece of code disables updatin at time=0.0
%                 hasM = 0;
%                 return;
%             end
            hasM = (abs(round(this.time / this.measurementFreq) - (this.time / this.measurementFreq)) < sqrt(eps));
        end
        
        function [m] = measure(this)
            h = this.measureOp();
            m = h(this.truth) + this.measureErr();
        end
        
        function [h] = measureOp(this)
            h = this.h;
        end
        
        function [merr] = measureErr(this, N, rndStr)
            if nargin<2
                N = 1;
            end
            % we simulate measurement of the truth, so we take our "own" PRNG
            if nargin<3
                rndStr = this.evidenceNoiseRng;
            end
            merr = (chol(this.R) * randn(rndStr, this.measurementDimension, N));
        end
        
        function [cov] = measureCov(this)
            cov = this.R;
        end
        
        function [L] = distanceMatrix(this)
            error('RandomWalk does not work with localisation');
        end
        
        function [tCount] = timeStepCount(this)
            tCount = this.tCount;
        end
        
        [this] = plot(this, ens)
    end
end

