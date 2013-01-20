classdef Advection < models.Model
    %ADVECTION Simple linear advection model, see [Burgers1998] etc.
    
    properties (SetAccess = private)
        deterministicDimension;
        stochasticDimension;
        measurementDimension;
        decorrelationLength;
        measurementFreq;
        tCount;
        tStep;
        measurementStdDev;
        evidenceNoiseRng;
        h;
    end
    
    properties (Access = private)
        time = 0;
        H;
        R;
        L;
    end
    
    methods (Access = private)
        function H = computeH(this)
            H = spalloc(this.measurementDimension,this.deterministicDimension,this.measurementDimension);
            
            for i = 1:this.measurementDimension
                H(i, mod(floor(this.deterministicDimension/this.measurementDimension/2) + floor(this.deterministicDimension/this.measurementDimension)*i, this.deterministicDimension)) = 1; %#ok<SPRIX>
            end
        end
        
        function R = computeR(this)
            R = speye(this.measurementDimension) * this.measurementStdDev.^2;
        end
        
        function L = computeL(this)
            L = zeros(this.deterministicDimension, this.measurementDimension);
            
            % we use the positions where the measurement operator has its
            % maximum as measurement locations - note that this only works
            % if that approximation makes any sense! Otherwise consider
            % methods like described in Fertig2007a, appendix A&B.
            mpos = tools.computeMPos(this.measurementDimension, this.deterministicDimension, this.h);
            for i = 1:this.deterministicDimension
                for j = 1:this.measurementDimension
                    L(i,j) = min([norm(i - mpos(j)), norm(i - mpos(j) - this.deterministicDimension), norm(i - mpos(j) + this.deterministicDimension)]);
                end
            end
        end
    end
    
    methods (Access = protected)
        function state = createTruth(this)
            state = zeros(this.deterministicDimension,1);
            s = (1:this.deterministicDimension)';
            
            for i = 0:this.stochasticDimension/2-1
                wavenumber = i * 2.0 * pi / this.deterministicDimension;
                amplitude = rand;
                phaseshift = rand * 2.0 * pi;
                state = state + amplitude * sin(wavenumber*s + phaseshift);
            end
            state = state + 6; % offset
        end
        
        function x = defaultMeasurementOperator(this, representation)
            if (isa(representation, 'numeric')) % this is the case where we measure the "truth"
                x = this.H*representation;
            elseif (isa(representation, 'representations.Ensemble'))
                x = this.H*representation.ensembleM;
%             elseif (isa(representation, 'representations.PCE'))
%                 x = this.H*representation.coefficients;
            else
                error('defaultMeasurementOperator(): unsupported representation');
            end
        end
    end
    
    methods
        function this = Advection(varargin)
            p = inputParser;
            
            p.addParamValue('deterministicDimension', 1000, @(x)x>1);
            p.addParamValue('decorrelationLength', 25.0, @(x)x>1);
            p.addParamValue('measurementDimension', 4, @(x)x>=1);
            p.addParamValue('measurementFreq', 5, @(x)x>0);
            p.addParamValue('tCount', 300, @(x)x>0);
            p.addParamValue('tStep', 1, @(x)x>0);
            p.addParamValue('measurementStdDev', 0.1, @(x)x>0);
            p.addParamValue('measurementOperator', @(x) this.defaultMeasurementOperator(x), @(f) isa(f,'function_handle'));
            p.addParamValue('evidenceNoiseRng', RandStream('mlfg6331_64'), @(x) isa(x,'RandStream'));
            
            p.parse(varargin{:});
            
            this.deterministicDimension = p.Results.deterministicDimension;
            this.decorrelationLength = p.Results.decorrelationLength;
            this.measurementDimension = p.Results.measurementDimension;
            this.measurementFreq = p.Results.measurementFreq;
            this.tCount = p.Results.tCount;
            this.tStep = p.Results.tStep;
            this.measurementStdDev = p.Results.measurementStdDev;
            this.evidenceNoiseRng = p.Results.evidenceNoiseRng;
            this.h = p.Results.measurementOperator;
            
            assert(this.measurementFreq >= this.tStep, 'The measurement frequency cannot be smaller than the time step!');
            assert(mod(this.measurementFreq, this.tStep) < eps, 'The measurement frequency must be divisible by the time step!');
            
            this.stochasticDimension = floor(this.deterministicDimension/(2*this.decorrelationLength))*2;
            this.H = this.computeH();
            this.R = this.computeR();
            this.L = this.computeL();
            
            % create truth case
            this.truth = this.createTruth();
        end
        
        function [coefficients] = createInitialPCE(this, basis) %#ok<INUSD,MANU,STOUT>
            error('Advection does not work with PCE (yet...)');
        end
        
        function [] = createInitialEnsemble(this, ens, varargin)
%         function [state] = createInstances(this, sourceSamples)
            assert(isa(ens, 'representations.Ensemble'), 'createInitialEnsemble(): you must pass an Ensemble to this function');
            assert(size(ens.sampleM, 1) == this.stochasticDimension, 'createInitialEnsemble(): wrong stochastic dimension')

            N = size(ens.sampleM, 2);
            
            firstguess = this.truth + this.createTruth() - 6; % subtract offset
            
            state = zeros(this.deterministicDimension,N);
            s = repmat((1:this.deterministicDimension)',1,N);
            
            maxwave = this.stochasticDimension/2-1;
            
            pos = 1;
            
            for i = 0:maxwave
                wavenumber = i * 2.0 * pi / this.deterministicDimension;
                amplitude = repmat(ens.sampleM(pos,:), this.deterministicDimension,1);
                phaseshift = repmat(ens.sampleM(pos+1,:) * 2.0 * pi, this.deterministicDimension,1);
                state(:,:) = state(:,:) + amplitude .* sin(wavenumber*s + phaseshift);
                pos = pos+2;
            end
            
            state = bsxfun(@plus, state, firstguess-6);
            
            % rescale initial ensemble for exact variance of
            % 1.0 in each gridpoint
            m = mean(state, 2);
            
            state = bsxfun(@minus, state, m);
            state = bsxfun(@rdivide, state, std(state, 0, 2));
            state = bsxfun(@plus, state, m);
            state = bsxfun(@plus, state, 6); % offset
            ens.ensembleM = state;
        end
        
        function this = step(this, representation, t)
            assert(isa(representation, 'representations.Ensemble'), 'representation must be Ensemble');
            if nargin<2 || isempty(t)
                t = this.tStep;
            end
            assert(t >= 1);
            this.truth = this.truth( mod((1:end)-2, end)+1 );
            representation.ensembleM = representation.ensembleM( mod((1:end)-2, end)+1, :);
            this.time = this.time + 1;
        end
        
        function [err] = stepErr(this, t)
            assert(false, 'not implemented');
        end
        
        function [hasM] = hasMeasurement(this)
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
            L = this.L;
        end
        
        function [tCount] = timeStepCount(this)
            tCount = this.tCount;
        end
        
        [this] = plot(this, ens)
    end
end

