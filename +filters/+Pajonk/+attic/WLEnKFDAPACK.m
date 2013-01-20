classdef WLEnKFDAPACK < filters.Filter
    %WLENKF Wavelet-space localized ensemble Kalman filter (WLEnKF)
    
    properties
        sampling; % Type of sampling used for the measurement ensemble.
        sampleBias; % Either 'remove' or 'keep', defaults to 'remove'.
        samplingFactor; % The advanced sampling factor
        waveletLevel; % The amount of transformations (e.g. 3)
        expansionBase; % The scaling factor for different scale levels, e.g. 1.05 for 5%
        waveletMode; % Padding mode for the DWT
        wavelet; % The type of wavelets which are used for transformation
        wLocalizationFunction; % Function handle to the wavelet space localization function.
        mLocalizationFunction; % Function handle to the measurement space localization function.
        denoise; % if true, then a denoising of the ensemble members is performed
        
        WL;
        IWL;
        L;
    end
    
    methods
        function this = WLEnKFDAPACK(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('sampling', 'default', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            
            p.addOptional('samplingFactor', 1.0, @(x)x>=1.0 );
            p.addOptional('waveletLevel', 3, @(x)x>=1 );
            p.addOptional('expansionBase', 1.0, @(x)x>=0.0 );
            p.addOptional('denoise', false, @(x) isa(x, 'logical'));
            p.addOptional('wavelet', 'haar', ...
                @(x)any(strcmpi(x,{...
                'haar','db1','db2','db3','db4','db5','db6','db7','db8',...
                'db9','db10','db11','db12','db13','db14','db15','db16',...
                'db17','db18','db19','db20','db21','db22','db23','db24',...
                'db25','db26','db27','db28','db29','db30','db31','db32',...
                'db33','db34','db35','db36','db37','db38','db39','db40',...
                'db41','db42','db43','db44','db45',...
                'sym2','sym3','sym4','sym5','sym6','sym7','sym8',...
                'sym9','sym10','sym11','sym12','sym13','sym14','sym15','sym16',...
                'sym17','sym18','sym19','sym20','sym21','sym22','sym23','sym24',...
                'sym25','sym26','sym27','sym28','sym29','sym30','sym31','sym32',...
                'sym33','sym34','sym35','sym36','sym37','sym38','sym39','sym40',...
                'sym41','sym42','sym43','sym44','sym45',...
                'bior1.1', 'bior1.3', 'bior1.5',...
                'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8',...
                'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7',...
                'bior3.9', 'bior4.4', 'bior5.5', 'bior6.8',...
                'coif5',...
                })));
            p.addOptional('waveletMode', 'sym', ...
                @(x)any(strcmpi(x,{'zpd','sym','symh','asym','asymh','asymw','spd','sp1','ppd'})));
            p.addOptional('sampleBias', 'remove', ...
                @(x)any(strcmpi(x,{'keep','remove'})));
            
            p.addOptional('wLocalizationFunction',...
                @(L, dl, wl) exp(-1.0 .* (L./((1.0/wl)*1.0*dl)).^2),...
                @(f) isa(f, 'function_handle'));
            p.addOptional('mLocalizationFunction',...
                @(L, dl) exp(-1.0 .* (L./(1.0*dl)).^2),...
                @(f) isa(f, 'function_handle'));
            p.parse(varargin{:});
            
            this.sampling = p.Results.sampling;
            this.samplingFactor = p.Results.samplingFactor;
            this.sampleBias = p.Results.sampleBias;
            this.waveletLevel = p.Results.waveletLevel;
            this.expansionBase = p.Results.expansionBase;
            this.waveletMode = p.Results.waveletMode;
            this.wavelet = p.Results.wavelet;
            this.wLocalizationFunction = p.Results.wLocalizationFunction;
            this.mLocalizationFunction = p.Results.mLocalizationFunction;
            this.denoise = p.Results.denoise;
            
            dwtmode(this.waveletMode); % Mode how to pad signal. The
            % extension modes represent different ways of handling the
            % problem of border distortion in signal and
            % image analysis. WARNING: THIS SHOULD NEVER BE CHANGED DURING
            % COMPUTATION!!!
        end
        
        function this = update(this, model)
%             S = model.ensemble;
            [X, sizes, decompositionLength] = model.wavedec(model.ensemble, this.waveletLevel, this.wavelet);
            
            % compute default localization (used for measurements)
            if isempty(this.L)
                this.L = this.mLocalizationFunction(model.distanceMatrix(), model.decorrelationLength());
            end;
            
%             % compute wavelet space localization operator from wavelets and
%             % distances
%             if isempty(this.WL)
%                 this.WL = this.computeWaveletSpaceLocalizationOperatorEx(Lo_D,Hi_D,Lo_R,Hi_R,model.distanceMatrix(),model.decorrelationLength());
% %                 this.WL = this.computeWaveletSpaceLocalizationOperator2(Lo_D,Hi_D,Lo_R,Hi_R,model.distanceMatrix(),model.decorrelationLength(),lengths);
%                 % For debugging & visualization: the Wavelet space localization
%                 % operator can be transformed back, of course:
%                 this.IWL = zeros(size(this.L));
%                 for i = 1:size(this.WL,2)
%                    s = waverec(this.WL(:,i),lengths,Lo_R,Hi_R);
%                    this.IWL(:,i) = s;
%                 end
%             end
            
%             WLex = this.computeWaveletSpaceLocalizationOperatorEx(Lo_D,Hi_D,Lo_R,Hi_R,L,dL);
%             WL =
%             ;
            
            % Compute the measurement operator in Wavelet space. 
%             HW = computeWaveletSpaceMeasurementOperator(this,Lo_D,Hi_D,model.measureOp());
            % use the wavelet space measurement operator to compute the
            % ensemble of measurements:
%             HX = HW'*X;
            % "proof by example" for this fact:
            %test = HW'*X;
            %test2 = H'*S;
            %plot(test-test2);
            % BUT: As the Wavelet transform is linear, this can also be
            % omitted and the standard measurement operator of the model
            % can be used (see below). IMPORTANT: If the measurement
            % operator is non-linear, this does not work, of course!!
            
            % Compute simulated measurements directly from the
            % ensemble, as this is equivalent to the computation in Wavelet
            % space, but computationally easier.
            HX = model.ensMeasure();

            % prepare noisy measurement of the truth to assimilate
            d = model.measure();
            
            R = model.measureCov();
            H = model.measureOp();
            N = size(X, 2);
            
            % compute simulated measurement mean
            EHX = mean(HX,2);
            
            % compute ensemble mean
            EX = mean(X, 2);
            
            A = zeros(size(X));
            HA = zeros(size(HX));
            % substract effect mean from effects on ensemble
            for j=1:N
                A(:,j) = X(:,j) - EX;
                HA(:,j) = HX(:,j) - EHX;
            end
            
            % sample from measurement error
            if strcmp(this.sampling, 'advanced') && this.samplingFactor > 1.0
                % Advanced sampling
                E = zeros(length(d), this.samplingFactor*N);
                for j=1:this.samplingFactor*N
                    E(:,j) = model.measureErr();
                end
                
                E = tools.sampling.condenseSample(E, N);
            else
                % Default sampling
                E = zeros(length(d), N);
                for j=1:N
                    E(:,j) = model.measureErr();
                end
            end
            
            if strcmp(this.sampleBias, 'remove')
                % OPTIONAL: Remove sampling bias
                E = tools.sampling.removeBias(E);
            end
            
            % create measurement ensemble
            D = zeros(length(d), N);
            for j=1:N
                D(:,j) = E(:,j) + d;
            end
            
            % perform the Ensemble Kalman Filter update step
%             P = R + ((H*this.L).*((1/(N-1)) * (HA * HA')));
            P = R + (((1/(N-1)) * (HA * HA')));
            
            Pinv = inv(P);
            Z = HA' * Pinv; %#ok<MINV>
            
            K = (1 / (N-1)) * A * Z;
%             K = this.WL.*K; % localize K
            
%             wK = this.waveletLocalizedK(Lo_D,Hi_D,Lo_R,Hi_R,lengths,model.distanceMatrix(),model.decorrelationLength(),X,H,N,Z);
            
            % update
%             X = X + wK * (D - HX);
            X = X + K * (D - HX);
            
            % transform back
%             for i = 1:size(X,2)
%                 s = waverec(X(:,i),lengths,Lo_R,Hi_R);
%                 S(:,i) = s;
%             end
            model.ensSet(model.waverec(this.wavelet, X, sizes, decompositionLength));
        end
        
        % this algorithm is EVIL but works - but consider reimplementation
        function K = waveletLocalizedK(this, Lo_D,Hi_D,Lo_R,Hi_R,lengths, L, dL, X,H,N,Z)
            wX = zeros([lengths(end) size(X,2) this.waveletLevel+1]);
            wA = zeros([lengths(end) size(X,2) this.waveletLevel+1]);
            wK = zeros([lengths(end) size(H,1) this.waveletLevel+1]);
            Li = zeros([lengths(end) size(H,1) this.waveletLevel+1]);
            
            % compute a full reconstruction from every wavelet level
            for level = 1:this.waveletLevel+1
                for n = 1:N
                    wX(:,n,level) = wrcoef('a',X(:,n),lengths,Lo_R,Hi_R,level-1);
                end
            end
            
            % compute and subtract mean
            wEX = squeeze(mean(wX, 2));
            for n = 1:N
                wA(:,n,:) = squeeze(wX(:,n,:)) - wEX;
            end
            
            % compute kalman gain and localize per level
            % smoothly combine to multilevel gain cK
            cK = zeros(size(wK,1), size(wK,2));
            clcorr = zeros(size(wK,1), size(wK,2));
            for level = 1:this.waveletLevel+1
%                 Li = this.wLocalizationFunction(L, this.expansionBase^(level) * dL, this.waveletLevel-level+2);
                Li(:,:,level) = this.wLocalizationFunction(L, dL, this.waveletLevel-level+2);
                wK(:,:,level) = ((1 / (N-1)) * squeeze(wA(:,:,level)) * Z);
                
                if (level == 1)
                    cK = squeeze(Li(:,:,level)) .* wK(:,:,level);
                    clcorr = squeeze(Li(:,:,level)) .* ((1 / (N-1)) * squeeze(wA(:,:,level)) * (H*squeeze(wA(:,:,level)))');
                else
                    cK = cK + (squeeze(Li(:,:,level)) - squeeze(Li(:,:,level-1))) .* wK(:,:,level);
                    clcorr = clcorr + (squeeze(Li(:,:,level)) - squeeze(Li(:,:,level-1))) .* (((1 / (N-1)) * squeeze(wA(:,:,level)) * (H*squeeze(wA(:,:,level)))'));
                end
            end
            
            defCorr = Li(:,:,this.waveletLevel+1) .* ((1 / (N-1)) * squeeze(wA(:,:,1)) * (H*squeeze(wA(:,:,1)))');
            L = Li(:,:,this.waveletLevel+1);
             
%             figure(18);
%             plot(cK);
%             figure(19);
%             plot(clcorr);
%             figure(1);
            
            K = zeros(size(X,1),size(cK,2));
            for i = 1:size(K,2)
                K(:,i) = wavedec(cK(:,i), this.waveletLevel,Lo_D,Hi_D);
            end
            
           
        end
        
        % This function can be used to optionally transform the linear
        % measurement operator H to Wavelet space.
        function HW = computeWaveletSpaceMeasurementOperator(this,Lo_D,Hi_D,H)
            H = full(H)';
            [coefficients, ~] = wavedec(H(:,1),this.waveletLevel,Lo_D,Hi_D);
            HW = zeros(length(coefficients), size(H,2));
            for i = 1:size(H,2)
                [coefficients, ~] = wavedec(H(:,i),this.waveletLevel,Lo_D,Hi_D);
                HW(:,i) = coefficients;
            end
        end
        
        function locOp = computeWaveletSpaceLocalizationOperator(this,Lo_D,Hi_D,Lo_R,Hi_R,L,dL)
            measurements = size(L,2);
            
            locOp = [];
            
            for level = 1:this.waveletLevel
                % New variable for wavelet-space localization:
                % Localization function for different levels can be
                % different! First approach: same function, but double
                % the effective decorrelation length for coarser levels, as
                % correlations between slowly varying modes are expected to
                % "make physical sense" between more distant points in
                % space.
%                 Li = this.wLocalizationFunction(L, 2^(level-1) * dL);
                Li = this.wLocalizationFunction(L, dL, level);
                
                % Construction of localization operator of *correct length*
                % by Wavelet analysis. This is definitely not the most
                % efficient method, but the most obvious one.
                levelCoef = [];
                for measurement = 1:measurements
                    [coefficients,lengths] = wavedec(Li(:,measurement),this.waveletLevel,Lo_D,Hi_D);
                    coef = appcoef(coefficients,lengths,Lo_R,Hi_R,level);
                    % Normalize the filter coefficients to be at most 1.0.
                    % This is absolutely necessary to avoid "explosion".
                    coef = coef ./ max(coef);
                    
%                     if (level > 4)
%                         levelCoef = cat(2, levelCoef, zeros(size(coef)));
%                     else
                        levelCoef = cat(2, levelCoef, coef);
%                     end
                end
                
                locOp = cat(1, levelCoef, locOp);
                
                if (level == this.waveletLevel)
                    locOp = cat(1, levelCoef, locOp);
                end
            end
        end
        
        
        function locOp = computeWaveletSpaceLocalizationOperatorEx(this,Lo_D,Hi_D,Lo_R,Hi_R,L,dL,expansionBase)
            measurements = size(L,2);
            
            locOp = [];
            
            for level = 1:this.waveletLevel
                % New variable for wavelet-space localization:
                % Localization function for different levels can be
                % different! First approach: same function, but double
                % the effective decorrelation length for coarser levels, as
                % correlations between slowly varying modes are expected to
                % "make physical sense" between more distant points in
                % space.
%                 Li = this.wLocalizationFunction(L, 2^(level-1) * dL, level);
                Li = this.wLocalizationFunction(L, this.expansionBase^(level) * dL, level);
                
                Li = [Li;Li;Li;Li;Li]; %#ok<AGROW>
                
                % Construction of localization operator of *correct length*
                % by Wavelet analysis. This is definitely not the most
                % efficient method, but the most obvious one.
                levelCoef = [];
                for measurement = 1:measurements
                    [coefficients,lengths] = wavedec(Li(:,measurement),this.waveletLevel,Lo_D,Hi_D);
                    coef = appcoef(coefficients,lengths,Lo_R,Hi_R,level);
                    % Normalize the filter coefficients to be at most 1.0.
                    % This is absolutely necessary to avoid "explosion".
                    coef = coef ./ max(coef);
                    levelCoef = cat(2, levelCoef, coef);
                end
                
                locOp = cat(1, levelCoef, locOp);
                
                if (level == this.waveletLevel)
                    locOp = cat(1, levelCoef, locOp);
                end
            end
            
            [m,n] = size(L);
            ilocOP = zeros(m*5,n);
            for i = 1:size(locOp,2)
               s = waverec(locOp(:,i),lengths,Lo_R,Hi_R);
               ilocOP(:,i) = s;
            end
            
            ilocOP = ilocOP(2*m+1:3*m,:);
            
            s =  wavedec(ilocOP(:,1),this.waveletLevel,Lo_D,Hi_D);
            clear locOp;
            locOp = zeros(length(s), n);
            for i = 1:size(ilocOP,2)
               s =  wavedec(ilocOP(:,i),this.waveletLevel,Lo_D,Hi_D);
               locOp(:,i) = s;
            end
        end
        
        function locOp = computeWaveletSpaceLocalizationOperator2(this,Lo_D,Hi_D,Lo_R,Hi_R,L,dL,lengths)
            measurements = size(L,2);
            
%             filter = exp(-1.0 .* (abs(linspace(-2.5, 2.5, waveletFilterSize)).^2));
            filter = Lo_D;
            
            locOp = [];
            
            for level = 1:this.waveletLevel
                Li = this.wLocalizationFunction(L, this.expansionBase^(level) * dL, level);
                
                Linew = [];
                for measurement = 1:measurements
                    m = Li(:,measurement);
                    for j = 1:level
                        m = dyaddown(conv(m, filter));
                    end
                    
                    Linew = cat(2, Linew, m);
                end
                
                Linew = Linew ./ (max(max(Linew)));
                
                locOp = cat(1, Linew, locOp);
                
                if (level == this.waveletLevel)
                    locOp = cat(1, Linew, locOp);
                end
            end
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

