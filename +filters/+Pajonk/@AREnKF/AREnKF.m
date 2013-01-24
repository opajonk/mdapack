classdef AREnKF < filters.Filter
    %ARENKF Advanced Regularisation ensemble Kalman filter (AREnKF)
    % In this EnKF, many modern regularisation methods are implemented and
    % can be tested mutually or in combination.
    
    properties
        opts; % structure for the options
        
        measurementEnsembleRng; % the PRNG used to generate the measurement ensemble
        
        N; % store ensemble size
        n; % store state space dimension
        m; % store measurement space dimension
    end
    
    
    properties(Access=private)
        Lo_D; % low-pass decomposition filter
        Hi_D; % high-pass decomposition filter
        Lo_R; % low-pass reconstruction filter
        Hi_R; % high-pass reconstruction filter
        wlX; % decomposition lengths for the ensemble (needed for reconstruction)
        wlHX; % decomposition lengths for the measurements (needed for reconstruction)
        
        % cache the localizers, since it may be expensive to compute them
        waveletGainLoc; 
        covLoc; 
        mCovLoc;
    end
    
    methods
        function this = AREnKF(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            %% solving
            p.addOptional('solving', 'eigen', ...
                @(x)any(strcmpi(x,{'svd','MATLAB','eigen'})));
            
            %% sampling
            p.addOptional('sampling_mode', 'standard', ...
                @(x)any(strcmpi(x,{'standard','advanced'})));
            p.addOptional('sampling_advanced_factor', 1.0, @(x)x>=1.0);
            p.addOptional('sampling_order', 2, @(x) any(find([0,1,2] == x)));
            
            %% screening
            p.addOptional('screening', 'none', ...
                @(x)any(strcmpi(x,{'none','gain','cov'})));
            p.addOptional('screening_function_gain',...
                'screening_gain_adaptivePajonkDetResA',...
                @(f) isa(f, 'char'));
            p.addOptional('screening_function_cov',...
                'screening_cov_adaptivePajonkDetResA',...
                @(f) isa(f, 'char'));
            
            %% localisation
            p.addOptional('localization', 'none', ...
                @(x)any(strcmpi(x,{'none','gain','measurement','grid','both'})));
            p.addOptional('localization_function',...
                @(L, dl, N) localizers.Gaspari1999_4_10(L, dl*1.3),...
                @(f) isa(f, 'function_handle'));
            % The GC4.10 is multiplied by 1.3 so that the width matches
            % approximately a Gaussian with stddev=dl.
            % Another possible localizer is the Gaussian: @(L, dl, N) exp(-1.0 .* (L./dl).^2);
            p.addOptional('localization_function_width',...
                1.0,...
                @(x)x>0.0);
            
            
            %% inflation
            p.addOptional('inflation', 'none', ...
                @(x)any(strcmpi(x,{'none','fixed','adaptive'})));
            p.addOptional('inflation_fixed_factor', 1.01, @(x)x>=1.0);
            
            %% wavelets
            p.addOptional('transformation_wavelet_enabled', false, @(x) isa(x, 'logical'));
            p.addOptional('transformation_wavelet_mode', 'grid', ...
                @(x)any(strcmpi(x,{'gain','both','measurement','grid'})));
            p.addOptional('transformation_wavelet_type', 'haar', ...
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
                'coif1','coif2','coif3','coif4','coif5',...
                })));
            p.addOptional('transformation_wavelet_level', Inf, @(x)x>=1 );
            p.addOptional('transformation_wavelet_extensionMode', 'per', ...
                @(x)any(strcmpi(x,{'zpd','sym','symh','asym','asymh','asymw','spd','sp1','ppd','per'})));
            p.addOptional('transformation_wavelet_localization_enabled', false, @(x) isa(x, 'logical'));
            p.addOptional('transformation_wavelet_localization_threshold_enabled', false, @(x) isa(x, 'logical'));
            p.addOptional('transformation_wavelet_localization_threshold', 0.2, @(x)x>=0.0);
            
            %% gridPoint normalization
            p.addOptional('normalization_gridPoint_enabled', false, @(x) isa(x, 'logical'));
            
            %% fourier normalization
            p.addOptional('normalization_fourier_enabled', false, @(x) isa(x, 'logical'));
            
            %% PRNGs
            p.addOptional('measurementEnsembleRng', ...
                RandStream('mt19937ar','Seed', now()), ...
                @(x) isa(x,'RandStream'));
            
            
            %% option parsing
            p.parse(varargin{:});
            this.opts = p.Results;
            
            % extract the PRNG from the options - it is not an option, so
            % remove it from there
            this.measurementEnsembleRng = this.opts.measurementEnsembleRng;
            this.opts = rmfield(this.opts, 'measurementEnsembleRng');
            
            if (this.opts.transformation_wavelet_enabled)
                % Mode how to pad signal. The
                % extension modes represent different ways of handling the
                % problem of border distortion in signal and
                % image analysis. WARNING: THIS SHOULD NEVER BE CHANGED DURING
                % COMPUTATION!!!
                dwtmode(this.opts.transformation_wavelet_extensionMode);
            end
            
            %% option sanity checks
            if (this.opts.transformation_wavelet_enabled)
                assert(strcmp(this.opts.localization,'none') || ~strcmp(this.opts.localization_function,'localization_gc410'), ...
                    'MDAPACK:invalidArgument',...
                    'localization_gc410 cannot be used with waveletTransform'...
                    );
            end
            if (this.opts.transformation_wavelet_enabled)
                assert(strcmp(this.opts.localization,'none') || ~strcmp(this.opts.localization_function,'localization_gaussian'), ...
                    'MDAPACK:invalidArgument',...
                    'localization_gaussian cannot be used with waveletTransform'...
                    );
            end
            if (this.opts.transformation_wavelet_localization_enabled)
                assert(this.opts.transformation_wavelet_enabled, ...
                    'MDAPACK:invalidArgument',...
                    'transformation_wavelet_localization_enabled must be used with transformation_wavelet_enabled'...
                    );
                assert(~strcmp(this.opts.solving,'MATLAB'), ...
                    'MDAPACK:invalidArgument',...
                    'transformation_wavelet_localisation_enabled cannot be used with the MATLAB solution method'...
                    );
            end
            
            assert(~strcmp(this.opts.localization,'gain'), ...
                'MDAPACK:invalidArgument',...
                'gain localisation not yet implemented'...
                );
            assert(~strcmp(this.opts.inflation,'adaptive'), ...
                'MDAPACK:invalidArgument',...
                'adaptive inflation not yet implemented'...
                );
            if (~strcmp(this.opts.screening,'none'))
                assert(~strcmp(this.opts.solving,'MATLAB'), ...
                    'MDAPACK:invalidArgument',...
                    'MATLAB solve cannot be used with Kalman gain screening - use svd or eigen'...
                    );
            end
            if (this.opts.normalization_fourier_enabled)
                warning('MDAPACK:invalidArgument', 'CAUTION: Fourier-space normalisation (white noise normalisation) is quite new and not fully implemented: the normalisation of the measurement noise is not fully implemented');
            end
            
            
            %% intialisations
            [this.Lo_D,this.Hi_D,this.Lo_R,this.Hi_R] = wfilters(this.opts.transformation_wavelet_type);
        end
        
        %% Update function (main functionality)
        this = update(this, model, representation)
        
        X = inflation(this, X)
        
        function str = char(this)
            str = class(this);
        end
    end
end

