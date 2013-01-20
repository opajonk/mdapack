classdef Lorenz1963 < models.Model
    %LORENZ1963 Deterministic non-periodic flow [Lorenz1963]
    
    properties (SetAccess = private)
        measurementSchedule;
        measurementDimension;
        tEnd;
        tStep;
        measurementStdDev;
        h;
        R;
        deterministicDimension;
        stochasticDimension;
    end
    
    properties (Access = private)
        time = 0;
        ensMeanData;
        ensStdData;
        summaryData
        truthData;
        mPosData;
        tUnit;
        timeScale;
        plotMode;
        plot2dDim;
        plotPDF;
        H;
        evidenceNoiseRng;
        firstGuess;
        firstGuessStdDev;
        firstGuessLogNormal;
        initialValues;
        measurementOperator;
        paper_pdftruthplot;
        ode45options;
    end
    
    methods %(Static = true)
        function [ dx ] = deterministicForward(this, ~, X)
            % Implementation of the ODEs of [Lorenz1963]
            X = reshape(X, this.stochasticDimension, []);
            
            if (this.stochasticDimension == 6)
                rho = X(4,:);
                sigma = X(5,:);
                beta = X(6,:);
            else
                rho = this.initialValues(4,1);
                sigma = this.initialValues(5,1);
                beta = this.initialValues(6,1);
            end
            
            dx = zeros(size(X));
            dx(1,:) = sigma .* (X(2,:) - X(1,:));
            dx(2,:) = X(1,:).*(rho - X(3,:)) - X(2,:);
            dx(3,:) = X(1,:).*X(2,:) - beta.*X(3,:);
            % dx(4:6,:) = 0; % parameters are constant, needs no update
            dx = reshape(dx, [], 1); % create column vector from solution
        end
    end
    
    methods %(Static = true)
        function [ dx ] = stochasticForward(this, ~, x, pce )
            % Implementation of the ODEs of [Lorenz1963]
            x = reshape(x, this.stochasticDimension, []);
            
            dx = zeros(size(x));
            
            if (this.stochasticDimension == 6)
                pcerho = x(4,:);
                pcesigma = x(5,:);
                pcebeta = x(6,:);
            else
                pcerho = zeros(1, size(pce.basis, 1));
                pcerho(1,1) = this.initialValues(4,1);
                
                pcesigma = zeros(1, size(pce.basis, 1));
                pcesigma(1,1) = this.initialValues(5,1);
                
                pcebeta = zeros(1, size(pce.basis, 1));
                pcebeta(1,1) = this.initialValues(6,1);
            end
            
            dx(1,:) = pce_multiply(pcesigma, pce.basis, (x(2,:) - x(1,:)), pce.basis, pce.basis, 'M', pce.multTensor);
            dx(2,:) = pce_multiply(x(1,:), pce.basis, (pcerho - x(3,:)), pce.basis, pce.basis, 'M', pce.multTensor) - x(2,:);
            dx(3,:) = pce_multiply(x(1,:), pce.basis, x(2,:), pce.basis, pce.basis, 'M', pce.multTensor) - pce_multiply(pcebeta, pce.basis, x(3,:), pce.basis, pce.basis, 'M', pce.multTensor);
            % dx(4:6,:) = 0; % parameters are constant, needs no update
            dx = reshape(dx, [], 1); % create column vector from solution
        end
    end
    
    methods (Access = private)
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
        
        function R = computeR(this)
            R = speye(this.measurementDimension) * this.measurementStdDev.^2;
        end
    end
    
    methods
        function this = Lorenz1963(varargin)
            p = inputParser;
            
            % !!!!!!!!!!!!!!!! IF tUnit IS CHANGED tStep IS INFLUENCED !!!!!!!!!!!!!!!!!!!!
            % !! IT MAY BE NECESSARY TO CHANGE IT TO REMAIN STABLE WITH THE ODE SOLUTION !!
            p.addParamValue('tUnit', 5, @(x)x>0); % time unit is 5 days, like in [Lorenz1984], [Lorenz1996]
            % However, [Lorenz1963] did not specify any time unit [Lorenz1963, p. 135, 136]
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            p.addParamValue('firstGuess', [3; -3.0; 20.0; 28; 10; 8/3]); % the first guess mean for the estimator (MUST match the stochastic dimension!)
            p.addParamValue('firstGuessStdDev', [1.0; 1.0; 1.0; 0; 0; 0]); % the std dev of the first guess mean for the estimator, 0.0 == not uncertain
            p.addParamValue('firstGuessLogNormal', false, @(x) isa(x,'logical')); % if yes, then the parameters have a lognormal distribution instead of a normal one
            % rho = 28; sigma = 10; beta = 8/3;
            % usual values from [Lorenz1963, p. 136]
            p.addParamValue('truth', [1.508870; -1.531271; 25.46091; 28; 10; 8/3], @(x) isa(x, 'double')); % Taken from [Evensen2009, p. 84]
            p.addParamValue('measurementSchedule', 1.25:1.25:100, @(x) isa(x,'double')); % measurement every 1.25 days
            p.addParamValue('measurementFreq', 1.25, @(x)x>0); % measurement interval = 0.25 ( * tUnit = 5!) taken from [Evensen2009, p. 87]
            p.addParamValue('tEnd', 100, @(x)x>0); % t=[0,20]*tUnit = [0,100] [Evensen2009, p. 87]
            p.addParamValue('tStep', 0.01, @(x)x>0); % tStep is chosen to make RK4 stable and create nice smooth plots ;-)
            p.addParamValue('measurementStdDev', sqrt(0.5), @(x)x>0); % variance = 0.5 from [Evensen2009, p. 84]
            p.addParamValue('plotMode', '2d', ...
                @(x)any(strcmpi(x,{'2d','3d'})));
            p.addParamValue('plot2dDim', 1, @(x) any(find([1,2,3,4,5,6] == x)));
            p.addParamValue('plotPDF', false, @(x) isa(x,'logical')); % if yes, plot a PDF estimate in 2D mode
            p.addParamValue('measurementOperator', @(x) this.defaultMeasurementOperator(x), @(f) isa(f,'function_handle'));
            p.addParamValue('evidenceNoiseRng', RandStream('mt19937ar','Seed',tools.shuffleSeed), @(x) isa(x,'RandStream'));
            
            % If yes, plot the truth and the measurement into a global figure
            % This is meant to be used with
            % controllers.singleRun_PDFPlot().
            p.addParamValue('paper_pdftruthplot', false, @(x) isa(x,'logical'));
            
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
            this.firstGuess = p.Results.firstGuess;
            this.firstGuessStdDev = p.Results.firstGuessStdDev;
            this.firstGuessLogNormal = p.Results.firstGuessLogNormal;
            this.truth = p.Results.truth;
            this.initialValues = p.Results.truth;
            this.measurementOperator = p.Results.measurementOperator;
            this.paper_pdftruthplot = p.Results.paper_pdftruthplot;
            
            
            this.ode45options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
            
            % Taken from [Lorenz1963, p. 137]
            % this.truth = [0;1;0];
            
            if (p.Results.firstGuessStdDev(4,1) > 0.0)
                % we have a 6-D system
                this.truth = this.initialValues(1:6,1);
                
                % in case the default measurement operator is used we need this H-matrix
                % do not measure the parameters - this is not possible
                this.H = [1 0 0 0 0 0; ...
                          0 1 0 0 0 0; ...
                          0 0 1 0 0 0];
                
                this.deterministicDimension = 6;
                this.stochasticDimension = 6;
                
                this.firstGuess = p.Results.firstGuess(1:6,1);
                this.firstGuessStdDev = p.Results.firstGuessStdDev(1:6,1);
                
            else
                % we have a 3-D system
                this.truth = this.initialValues(1:3,1);
                
                % in case the default measurement operator is used we need this H-matrix
                this.H = [1 0 0; ...
                          0 1 0; ...
                          0 0 1];
                this.deterministicDimension = 3;
                this.stochasticDimension = 3;
                
                this.firstGuess = p.Results.firstGuess(1:3,1);
                this.firstGuessStdDev = p.Results.firstGuessStdDev(1:3,1);
                
            end
            
            this.measurementDimension = length(this.h(this.truth)); % compute the dimension of the measurement
            this.R = this.computeR();
            
%             assert(max(this.measurementSchedule) <= this.tEnd, 'The measurement schedule cannot last longer than tEnd');
            assert(min(this.measurementSchedule) >= 0, 'The measurement schedule cannot start before 0');
            assert(isempty(find(arrayfun(@(x) mod(x, this.tStep) > eps, this.measurementSchedule), 1)), 'All entries in the schedule must be divisible by the time step');
        end
        
        function [] = createInitialPCE(this, pce)
            assert(isa(pce, 'representations.PCE'), 'you must pass a PCE to this function');
            %             firstguess = this.measureErr();
            
            if (this.stochasticDimension == 6)
                firstguess = this.firstGuess;
                pce.coefficients = zeros(size(pce.basis')); % PCE coefficients must be row vector
                
                if this.firstGuessLogNormal
                    pce.coefficients(1:3,1) = firstguess(1:3); % mean
                    pce.coefficients(1:3,2:4) = diag(this.firstGuessStdDev(1:3)); % first order terms, all others zero
                    
                    % The parameters are lognormal. Attention: the firstguess and
                    % firstGuessStdDev are the parameters for the *underlying normal
                    % distribution*, not for the lognormal one!
                    % A formula for the PCE coefficients of lognormal RV
                    % is given by Ghanem1999, eq. (45).
                    
                    for i = 4:6
                        idx = 1:6;
                        idx(i) = [];
                        
                        % find those PCE basis positions, where all
                        % non-i-th RV components are zero.
                        % These are the coefficients which describe
                        % independent priors, and therefore we need to set
                        % those.
                        pos = sum(pce.basis(:,idx),2) == 0; 
                        pce.coefficients(i,pos) = exp(firstguess(i) + this.firstGuessStdDev(i)^2/2) * this.firstGuessStdDev(i).^(pce.basis(pos,i))./pce.multFactorial(pos)';
                    end
                else
                    pce.coefficients(:,1) = firstguess; % mean
                    pce.coefficients(1:6,2:7) = diag(this.firstGuessStdDev); % first order terms, all others zero
                end
            else
                firstguess = this.firstGuess(1:3,1);
                pce.coefficients = zeros(size(pce.basis')); % PCE coefficients must be row vector
                pce.coefficients(:,1) = firstguess; % mean
                pce.coefficients(:,2:4) = diag(this.firstGuessStdDev); % first order terms, all others zero
            end
            
            assert(isequal(size(pce.coefficients), size(pce.basis')),'createInitialPCE(): dimensions of basis and coefficients do not match');
        end
        
        function [] = createInitialEnsemble(this, ens)
            assert(isa(ens, 'representations.Ensemble'), 'createInitialEnsemble(): you must pass an Ensemble to this function');
            % sourceSamples should be standard Gaussian distributed
            assert(size(ens.sampleM, 1) == this.stochasticDimension, 'wrong stochastic dimension')
            
            % we create the instances using PC to be able to compare PC forward and MC forward
            if (this.stochasticDimension == 6)
                firstguess = this.firstGuess(1:6,1);
                
                if this.firstGuessLogNormal
                    % use PCE of order 5 to have a reasonable approximation of the lognormal RV
                    basis = multiindex(this.stochasticDimension, 6); 
                    multFactorial = multiindex_factorial(basis)';
                    coefficients = zeros(size(basis')); % PCE coefficients must be row vector
                    
                    coefficients(1:3,1) = firstguess(1:3); % mean
                    coefficients(1:3,2:4) = diag(this.firstGuessStdDev(1:3)); % first order terms, all others zero
                    
                    % The parameters are lognormal. Attention: the firstguess and
                    % firstGuessStdDev are the parameters for the *underlying normal
                    % distribution*, not for the lognormal one!
                    % A formula for the PCE coefficients of lognormal RV
                    % is given by Ghanem1999, eq. (45).
                    
                    for i = 4:6
                        idx = 1:6;
                        idx(i) = [];
                        
                        % find those PCE basis positions, where all
                        % non-i-th RV components are zero.
                        % These are the coefficients which describe
                        % independent priors, and therefore we need to set
                        % those.
                        pos = sum(basis(:,idx),2) == 0; 
                        coefficients(i,pos) = exp(firstguess(i) + this.firstGuessStdDev(i)^2/2) * this.firstGuessStdDev(i).^(basis(pos,i))./multFactorial(pos)';
                    end
                else
                    basis = multiindex(this.stochasticDimension, 1); % PCE of order 1 is enough
                    coefficients(:,1) = firstguess; % mean
                    coefficients(1:6,2:7) = diag(this.firstGuessStdDev); % first order terms, all others zero
                end
            else
                firstguess = this.firstGuess(1:3,1);
                basis = multiindex(this.stochasticDimension, 1); % PCE of order 1 is enough
                
                coefficients = zeros(size(basis')); % PCE coefficients must be row vector
                coefficients(:,1) = firstguess; % mean
                coefficients(1:3,2:4) = diag(this.firstGuessStdDev(1:3,1)); % first order terms, all others zero, like [Shen2010]
            end
            
            ens.ensembleM = pce_evaluate(coefficients, basis, ens.sampleM);
        end
        
        function this = step(this, representation, t)
            if nargin<2 || isempty(t)
                t = this.tStep;
            end
            
            endTime = this.tStep / this.tUnit;
            
            [~,Y] = ode113(@(t, X) this.deterministicForward(t, X), [0 endTime], this.truth', this.ode45options);
            this.truth = Y(end,:)';
                        
            % advance the representation - different implementations may be possible here
            if isa(representation, 'representations.Ensemble')
                [~,Y] = ode113(@(t, X) this.deterministicForward(t, X), [0 endTime], representation.ensembleM, this.ode45options);
                representation.ensembleM = reshape(Y(end,:), size(representation.ensembleM));
            elseif isa(representation, 'representations.PCE') || isa(representation, 'representations.KLEPCE')
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
            
            if(this.paper_pdftruthplot)
                global ha;
                % global fig;
                
                plot(ha, [this.truth(1)],[0],'-kx','LineWidth', 2.0, 'MarkerSize', 15);
                % plot(ha, [tru(1) tru(1)],[0 1000],'-bo','LineWidth', 1.0);
                % plot(ha(2), [tru(2) tru(2)],[0 1000],'-b','LineWidth', 1.0);
                % plot(ha(3), [tru(3) tru(3)],[0 1000],'-b','LineWidth', 1.0);
                
%                 plot(ha, [m(1)],[0],'-b+','LineWidth', 2.0, 'MarkerSize', 15);
                % plot(ha(2), [m(2) m(2)],[0 1000],'-g','LineWidth', 1.0);
                % plot(ha(3), [m(3) m(3)],[0 1000],'-g','LineWidth', 1.0);
            end
        end
        
        function [f] = measureOp(this)
            f = this.measurementOperator;
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
            assert(false, 'Lorenz1963 does not support localization');
        end
        
        function [dl] = decorrelationLength(this)
            assert(false, 'Lorenz1963 does not support localization');
        end
        
        function [tCount] = timeStepCount(this)
            tCount = this.tEnd / this.tStep;
        end
        
        [this] = plot(this, representation, varargin)
    end
end


