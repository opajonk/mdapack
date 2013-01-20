classdef Lorenz1996 < models.Model
    % LORENZ1996 The spatial model of [Lorenz1996], [Lorenz2006]
    
    properties (SetAccess = private)
        measurementSchedule;
        tEnd;
        tStep;
        measurementStdDev;
        deterministicDimension;
        stochasticDimension
        measurementDimension;
        decorrelationLength;
        measurementResolution;
    end
    
    properties (Access = private)
        time = 0;
        ensMeanData;
        ensStdData;
        truthData;
        timeScale;
        H;
        R;
        L;
        forcing;
        mPosData;
        plotMode;
        plot2dDim;
        tUnit;
        evidenceNoiseRng;
        h;
        precomputedPrior;
        ode45options;
    end
    
    methods (Static = true)
        function [ dx ] = deterministicForward( ~, X, F)
            % Implementation of the ODEs of Lorenz2006 (Lorenz1996) with
            % values taken from Lorenz2006, p. 44
            ds = size(X,1);
            dx = zeros(size(X));
            
            % indices as helper arrays - cyclic domain
            km2 = [ds-1 ds 1:ds-2]; % k-2
            km1 = [ds 1:ds-1]; % k-1
            kp1 = [2:ds 1]; % k+1
            
            dx(:,:) = -1.0.*X(km2,:).*X(km1,:) + X(km1,:).*X(kp1,:) - X(:,:) + F;
        end
        
        function [ dx ] = stochasticForward( ~, coefficents, F, pce )
            x = reshape(coefficents, size(pce.basis'));
            
            pceF = zeros(1, size(pce.basis, 1));
            pceF(1,1) = F;
            
            ds = size(x,1);
            dx = zeros(size(x));
            
            
            dx(1,:) = -1.0 * pce_multiply(x(ds-1,:), pce.basis, x(ds,:), pce.basis, pce.basis, 'M', pce.multTensor) ...
                    + pce_multiply(x(ds,:), pce.basis, x(2,:), pce.basis, pce.basis, 'M', pce.multTensor) - x(1,:) + pceF;
                
            dx(2,:) = -1.0 * pce_multiply(x(ds,:), pce.basis, x(1,:), pce.basis, pce.basis, 'M', pce.multTensor) ...
                    + pce_multiply(x(1,:), pce.basis, x(3,:), pce.basis, pce.basis, 'M', pce.multTensor) - x(2,:) + pceF;
                
            for k = 3:ds-1
                dx(k,:) = -1.0 * pce_multiply(x(k-2,:), pce.basis, x(k-1,:), pce.basis, pce.basis, 'M', pce.multTensor) ...
                    + pce_multiply(x(k-1,:), pce.basis, x(k+1,:), pce.basis, pce.basis, 'M', pce.multTensor) - x(k,:) + pceF;
            end
            
            dx(ds,:) = -1.0 * pce_multiply(x(ds-2,:), pce.basis, x(ds-1,:), pce.basis, pce.basis, 'M', pce.multTensor) ...
                    + pce_multiply(x(ds-1,:), pce.basis, x(1,:), pce.basis, pce.basis, 'M', pce.multTensor) - x(ds,:) + pceF;
            
            dx = reshape(dx, [], 1); % create column vector from solution
        end
    end
    
    methods (Access = protected)
        function state = createTruth(this)
            ensFileName = sprintf('resources/Lorenz96.initEns.%d.mat', this.deterministicDimension);
            
            if (exist(ensFileName,'file'))
                % we can use the pre-computed states which are free of
                % transient effects
                s = load(ensFileName,'X');
                % the last ensemble member is the "truth"
                state = s.X(:,10000);
            else
                % we have to initialise randomly anyway
                % truth is ??? (Lorenz2006 says "rather arbitrary values of the variables"...)
                state = randn(this.deterministicDimension, 1);
            end
        end
    end
    
    methods (Access = private)
        function H = computeH(this)
            
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
        
        function x = defaultMeasurementOperator(this, representation)
            if (isa(representation, 'numeric')) % this is the case where we measure the "truth"
                x = this.H*representation;
            elseif (isa(representation, 'representations.Ensemble'))
                x = this.H*representation.ensembleM;
            elseif (isa(representation, 'representations.PCE'))
                x = this.H*representation.coefficients;
            else
                error('defaultMeasurementOperator(): unsupported representation');
            end
        end
    end
    
    methods
        function this = Lorenz1996(varargin)
            p = inputParser;
            
            % !!!!!!!!!!!!!!!!IF tUnit IS CHANGED tStep IS INFLUENCED !!!!!!!!!!!!!!!!!!!!!
            % !! IT MAY BE NECESSARY TO CHANGE IT TO REMAIN STABLE WITH THE ODE SOLUTION !!
            p.addParamValue('tUnit', 5, @(x)x>0); % time unit is 5 days (Lorenz2006, p. 44)
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            p.addParamValue('measurementSchedule', 2:2:365, @(x) isa(x,'double')); % measurement every 2 days
            p.addParamValue('measurementDimension', 20, @(x)x>=1); % measurements at every second grid point
            p.addParamValue('forcing', 8.0, @(x)x>=1); % the external forcing parameter. defaults to 8.0, as in [Lorenz2006, p. 44]
            p.addParamValue('decorrelationLength', 4.0, @(x)x>=1); % chosen on my own, based on the ODE
            p.addParamValue('measurementResolution', 0.0, @(x)x>=0); % 0.0 means Dirac operator (point measurement)
            p.addParamValue('deterministicDimension', 40, @(x)x>=4); % K > 3, [Lorenz2006, p. 44]
            p.addParamValue('tEnd', 365, @(x)x>0); % one year of integration
            p.addParamValue('tStep', 0.25, @(x)x>0); % time stepping = 0.05 time units = 6h = 0.25 days [Lorenz2006, p. 44]
            p.addParamValue('measurementStdDev', 1, @(x)x>0); % Accoring to Sakov2012 and others, this is a kind of "standard" setting
            p.addParamValue('measurementOperator', @(x) this.defaultMeasurementOperator(x), @(f) isa(f,'function_handle'));
            p.addParamValue('evidenceNoiseRng', RandStream('mt19937ar','Seed',tools.shuffleSeed), @(x) isa(x,'RandStream'));
            p.addParamValue('precomputedPrior', false, @(x) isa(x,'logical')); % if yes, a pre-computed ensemble is used for ensemble sizes <= 9999
            
            p.parse(varargin{:});
            
            this.measurementSchedule = p.Results.measurementSchedule;
            this.decorrelationLength = p.Results.decorrelationLength;
            this.measurementResolution = p.Results.measurementResolution;
            this.measurementDimension = p.Results.measurementDimension;
            this.forcing = p.Results.forcing;
            this.deterministicDimension = p.Results.deterministicDimension;
            this.stochasticDimension = p.Results.deterministicDimension;
            this.tEnd = p.Results.tEnd;
            this.tUnit = p.Results.tUnit;
            this.tStep = p.Results.tStep;
            this.measurementStdDev = p.Results.measurementStdDev;
            this.h = p.Results.measurementOperator;
            this.evidenceNoiseRng = p.Results.evidenceNoiseRng;
            this.precomputedPrior = p.Results.precomputedPrior;
            
            this.ode45options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
            
            % truth
            this.truth = this.createTruth();
            
            % in case the default measurement operator is used we need this H-matrix
            this.H = spalloc(this.measurementDimension,this.deterministicDimension,this.measurementDimension);
            for i = 0:this.measurementDimension-1
                this.H(i+1, round(mod((this.deterministicDimension/this.measurementDimension)*i, this.deterministicDimension)) + 1) = 1; 
            end
            
            this.measurementDimension = length(this.h(this.truth)); % compute the dimension of the measurement
            
            this.R = this.computeR();
            this.L = this.computeL();
            
            assert(max(this.measurementSchedule) <= this.tEnd, 'The measurement schedule cannot last longer than tEnd');
            assert(min(this.measurementSchedule) >= 0, 'The measurement schedule cannot start before 0');
            assert(isempty(find(arrayfun(@(x) mod(x, this.tStep) > eps, this.measurementSchedule), 1)), 'All entries in the schedule must be divisible by the time step');
        end
        
        function [] = createInitialPCE(this, pce)
            assert(isa(pce, 'representations.PCE'), 'createInitialPCE(): you must pass a PCE to this function');
            
            if (pce.pceOrder > 1)
                error('createInitialPCE(): running the Lorenz-96 model with a PCE order larger than one cannot be done at the moment');
            end
            
            firstguess = this.truth + randn(this.deterministicDimension, 1);
            pce.coefficients = zeros(size(pce.basis')); % PCE coefficients must be row vector
            pce.coefficients(:,1) = firstguess;
            pce.coefficients(1:this.stochasticDimension,2:this.stochasticDimension+1) = eye(this.stochasticDimension);
        end
        
        function [] = createInitialEnsemble(this, ens)
            assert(isa(ens, 'representations.Ensemble'), 'createInitialEnsemble(): you must pass an Ensemble to this function');
            
            ensFileName = sprintf('resources/Lorenz96.initEns.%d.mat', this.deterministicDimension);
            
            if (exist(ensFileName,'file') && ens.sampleSize < 10000 && this.precomputedPrior)
                s = load(ensFileName, 'X');
                % randomly choose from the ensemble WITHOUT replacement
                
                rp = ens.rndStream.randperm(9999);
                rp = rp(1:ens.sampleSize);
                
                ens.ensembleM = s.X(:,rp);
                ens.sampleM = [];
            else
                assert(size(ens.sampleM, 1) == this.stochasticDimension, 'createInitialEnsemble(): wrong stochastic dimension')
                % we create the instances using PC to be able to compare PC forward and MC forward
                basis = multiindex(this.stochasticDimension, 1);
                
                firstguess = this.truth + randn(this.deterministicDimension, 1);
                coefficients = zeros(size(basis')); % PCE coefficients must be row vector
                coefficients(:,1) = firstguess;
                coefficients(1:this.stochasticDimension,2:this.stochasticDimension+1) = eye(this.stochasticDimension);
                ens.ensembleM = pce_evaluate(coefficients, basis, ens.sampleM);
            end
        end
        
        function this = step(this, representation, t)
            if nargin<2 || isempty(t)
                t = this.tStep;
            end

            [~,Y] = ode113(@( t, X) models.Lorenz1996.deterministicForward(t, X, this.forcing), [0 this.tStep / this.tUnit], this.truth', this.ode45options);
            this.truth = Y(end,:)';
            
            endTime = this.tStep / this.tUnit;
            F = this.forcing;
            
            if isa(representation, 'representations.Ensemble')
                [~,Y] = ode113(@( t, X) models.Lorenz1996.deterministicForward(t, X, F), [0 endTime], representation.ensembleM, this.ode45options);
                representation.ensembleM = reshape(Y(end,:), size(representation.ensembleM));
            elseif isa(representation, 'representations.PCE')
                [~,Y] = ode113(@(t, X) models.Lorenz1996.stochasticForward(t, X, F, representation), [0 endTime], reshape(representation.coefficients, [], 1), this.ode45options);
                representation.coefficients = reshape(Y(end,:), size(representation.basis'));
            end
            
            this.time = this.time + this.tStep;
        end
        
        function [err] = stepErr(this, t)
            assert(false, 'not implemented');
        end
        
        function [hasM] = hasMeasurement(this)
            test = @(x) norm(x - this.time) < sqrt(eps);
            hasM = (any(arrayfun(test, this.measurementSchedule)));
        end
        
        function [m] = measure(this)
            f = this.measureOp();
            m = f(this.truth) + this.measureErr();
        end
        
        function [f] = measureOp(this)
%             f = @(x) this.defaultMeasurementOperator(x);
            f = this.h;
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
            tCount = this.tEnd / this.tStep;
        end
                
        [this] = plot(this, representation, varargin)
    end
end

