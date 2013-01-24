classdef Lorenz1984 < models.Model
    %LORENZ1984 'littlest' General Circulation Model [Lorenz1984]
    % 
    % I have verified this model agains Lorenz2005, plot (3) and it
    % produces precisely the same results. 2010-10-15, O.Pajonk
    
    properties (SetAccess = private)
        measurementSchedule;
        tEnd;
        tStep;
        measurementStdDev;
        h;
        R;
        measurementDimension;
        deterministicDimension = 3;
        stochasticDimension = 3;
    end
    
    properties (Access = private)
        time = 0;
        ensMeanData;
        ensStdData;
        truthData;
        summaryData;
        timeScale;
        tUnit;
        mPosData;
        plotMode;
        plot2dDim;
        plotPDF;
        H;
        multTensor;
        evidenceNoiseRng;
        firstGuessStdDev;
        lastMeasurement;
        trueForcing;
        firstGuess;
        initialValues;
        pceG;
        pceF;
        ode45options;
    end
    
    methods %(Static = true)
        function [ dx ] = deterministicForward(this, ~, X )
            X = reshape(X, this.stochasticDimension, []);
            % Implementation of the ODEs of Lorenz1984 with
            % values taken from Lorenz2005, p. 3f
            
            a = 0.25; b = 4.0;
            G = 1.23; F = 8.0;
            dx = zeros(size(X));
            dx(1,:) = -1.0*(X(2,:).^2) - (X(3,:).^2) - a*X(1,:) + a*F;
            dx(2,:) = X(1,:).*X(2,:) - b*X(1,:).*X(3,:) - X(2,:) + G;
            dx(3,:) = b*X(1,:).*X(2,:) + X(1,:).*X(3,:) - X(3,:);
            
            dx = reshape(dx, [], 1); % create column vector from solution
        end
        
        function [ dx ] = stochasticForward(this, ~, x, pce )
            x = reshape(x, this.stochasticDimension, []);
            % Implementation of the ODEs of Lorenz1984 with
            % values taken from [Lorenz2005], p. 3f and suggested by [Shen2010]
            a = 0.25; b = 4.0;
            dx = zeros(size(x));
            dx(1,:) = -1* pce_multiply(x(2,:), pce.basis, x(2,:), pce.basis, pce.basis, 'M', pce.multTensor) - pce_multiply(x(3,:), pce.basis, x(3,:), pce.basis, pce.basis, 'M', pce.multTensor) - a*x(1,:) + a*this.pceF;
            dx(2,:) = pce_multiply(x(1,:), pce.basis, x(2,:), pce.basis, pce.basis, 'M', pce.multTensor) - b*pce_multiply(x(1,:), pce.basis, x(3,:), pce.basis, pce.basis, 'M', pce.multTensor) - x(2,:) + this.pceG;
            dx(3,:) = b * pce_multiply(x(1,:), pce.basis, x(2,:), pce.basis, pce.basis, 'M', pce.multTensor) + pce_multiply(x(1,:), pce.basis, x(3,:), pce.basis, pce.basis, 'M', pce.multTensor) - x(3,:);
            dx = reshape(dx, [], 1); % create column vector from solution
        end
    end
    
    methods (Static = true)
        [ x ] = customMeasurementOperator(representation, H)
    end
    
    methods (Access = private)
        function R = computeR(this)
            % a diagonal (co-)variance matrix, containing the square of the
            % measurement standard deviation
            R = speye(this.measurementDimension) * this.measurementStdDev.^2;
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
        function this = Lorenz1984(varargin)
            p = inputParser;
            
            % !!!!!!!!!!!!!!!!IF tUnit IS CHANGED tStep IS INFLUENCED !!!!!!!!!!!!!!!!!!!!!
            % !! IT MAY BE NECESSARY TO CHANGE IT TO REMAIN STABLE WITH THE ODE SOLUTION !!
            p.addParamValue('tUnit', 5, @(x)x>0); % time unit is 5 days (Lorenz2005, p. 3)
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            p.addParamValue('firstGuess', [1.3; 0.3; -0.9]); % the first guess mean for the estimator (MUST match the stochastic dimension!)
            p.addParamValue('firstGuessStdDev', [0.5; 0.5; 0.5]); % if larger than 0, the forcing is considered as uncertain
            p.addParamValue('truth', [1.0; 0.0; -0.75], @(x) isa(x, 'double')); % truth is taken from Lorenz2005, p. 4 (truth for t=-300 --> his plot 3)
            p.addParamValue('tStep', 0.25, @(x)x>0); % time step = 0.25 days = 6h
            p.addParamValue('tEnd', 365, @(x)x>0); % integration for one year
            p.addParamValue('measurementSchedule', 2:2:365, @(x) isa(x, 'double')); % measurement every 2 days
            p.addParamValue('measurementStdDev', 0.1, @(x)x>=0);
            p.addParamValue('plotMode', '2d', ...
                @(x)any(strcmpi(x,{'2d','3d'})));
            p.addParamValue('plot2dDim', 1, @(x) any(find([1,2,3] == x)));
            p.addParamValue('plotPDF', false, @(x) isa(x,'logical')); % if yes, plot a PDF estimate in 2D mode
            p.addParamValue('measurementOperator', @(x) this.defaultMeasurementOperator(x), @(f) isa(f,'function_handle'));
            
            p.addParamValue('evidenceNoiseRng', RandStream('mt19937ar','Seed', now()), @(x) isa(x,'RandStream'));
            
            p.parse(varargin{:});
            
            this.measurementSchedule = p.Results.measurementSchedule;
            this.tEnd = p.Results.tEnd;
            this.tStep = p.Results.tStep;
            this.tUnit = p.Results.tUnit;
            this.measurementStdDev = p.Results.measurementStdDev;
            this.plotMode = p.Results.plotMode;
            this.plot2dDim = p.Results.plot2dDim;
            this.plotPDF = p.Results.plotPDF;
            this.h = p.Results.measurementOperator;            
            this.evidenceNoiseRng = p.Results.evidenceNoiseRng;
            this.initialValues = p.Results.truth;
            
            this.ode45options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
            
            
            this.truth = this.initialValues(1:3,1);
            
            % in case the default measurement operator is used we need this H-matrix
            this.H = [1 0 0 ; 0 1 0 ; 0 0 1 ];
            
            this.deterministicDimension = 3;
            this.stochasticDimension = 3;
            
            this.firstGuess = p.Results.firstGuess(1:3,1);
            this.firstGuessStdDev = p.Results.firstGuessStdDev(1:3,1);
            
            this.measurementDimension = length(this.h(this.truth)); % compute the dimension of the measurement
            this.R = this.computeR();
            
            assert(size(this.firstGuess, 1) == this.stochasticDimension, 'The size of the first guess and the stochastic dimension do not match');
%             assert(max(this.measurementSchedule) <= this.tEnd, 'The measurement schedule cannot last longer than tEnd');
            assert(min(this.measurementSchedule) >= 0, 'The measurement schedule cannot start before 0');
            assert(isempty(find(arrayfun(@(x) mod(x, this.tStep) > eps, this.measurementSchedule), 1)), 'All entries in the schedule must be divisible by the time step');
        end
        
        function [] = createInitialPCE(this, pce)
            assert(isa(pce, 'representations.PCE'), 'createInitialPCE(): you must pass a PCE to this function');
            assert(size(pce.basis, 2) == this.stochasticDimension, 'createInitialPCE(): wrong stochastic basis')
            
            this.pceG = zeros(1, size(pce.basis, 1));
            this.pceG(1,1) = 1.23;
            this.pceF = zeros(1, size(pce.basis, 1));
            this.pceF(1,1) = 8;
            
            firstguess = this.firstGuess(1:3,1);
            pce.coefficients = zeros(size(pce.basis')); % PCE coefficients must be row vector
            pce.coefficients(:,1) = firstguess; % mean
            pce.coefficients(1:3,2:4) = diag(this.firstGuessStdDev(1:3,1)); % first order terms, all others zero, like [Shen2010]
            
            assert(isequal(size(pce.coefficients), size(pce.basis')),'createInitialPCE(): dimensions of basis and coefficients do not match');
        end
        
        function [] = createInitialEnsemble(this, ens)
            assert(isa(ens, 'representations.Ensemble'), 'createInitialEnsemble(): you must pass an Ensemble to this function');
            assert(size(ens.sampleM, 1) == this.stochasticDimension, 'createInitialEnsemble(): wrong stochastic dimension')
            
            % we create the instances using PC to be able to compare PC forward and MC forward
            basis = multiindex(this.stochasticDimension, 1);
            
            firstguess = this.firstGuess(1:3,1);
            coefficients = zeros(size(basis')); % PCE coefficients must be row vector
            coefficients(:,1) = firstguess; % mean
            coefficients(1:3,2:4) = diag(this.firstGuessStdDev(1:3,1)); % first order terms, all others zero, like [Shen2010]
            ens.ensembleM = pce_evaluate(coefficients, basis, ens.sampleM);
        end
        
        function this = step(this, representation, t)
            if nargin<2 || isempty(t)
                t = this.tStep;
            end
            % create local copies of some variables to be able to use parfor
            endTime = this.tStep / this.tUnit;
            
            [~,Y] = ode113(@(t, X) this.deterministicForward(t, X), [0 endTime], this.truth', this.ode45options);
            this.truth = Y(end,:)';
            
            % advance the representation - different implementations may be possible here
            if isa(representation, 'representations.Ensemble')
                [~,Y] = ode113(@(t, X) this.deterministicForward(t, X), [0 endTime], representation.ensembleM, this.ode45options);
                representation.ensembleM = reshape(Y(end,:), size(representation.ensembleM));
            elseif isa(representation, 'representations.PCE')
                [~,Y] = ode113(@(t, X) this.stochasticForward(t, X, representation), [0 endTime], representation.coefficients, this.ode45options);
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
            
            if (norm(full(this.R)) > 0)
                merr = (chol(this.R) * randn(rndStr, this.measurementDimension, N));
            else
                merr = zeros(this.measurementDimension, N);
            end
        end
        
        function [cov] = measureCov(this)
            cov = this.R;
        end
        
        function [L] = distanceMatrix(this)
            assert(false, 'Lorenz1984 does not support localization');
        end
        
        function [dl] = decorrelationLength(this)
            assert(false, 'Lorenz1984 does not support localization');
        end
        
        function [tCount] = timeStepCount(this)
            tCount = this.tEnd / this.tStep;
        end
        
        [this] = plot(this, representation, varargin)
    end
end

