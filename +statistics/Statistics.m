classdef Statistics < handle
    %Statistics Default statistics class
    
    properties (Access = private)
        parms;
        model;
        results;
        filter;
        representation;
        ksdensityValues = 128;
    end
    
    methods
        function this = Statistics(parms, model, representation, filter)
            assert(isa(parms, 'struct'));
            assert(isa(model, 'models.Model'));
            assert(isa(representation, 'representations.Representation'));
            assert(isa(filter, 'filters.Filter'));
            
            this.parms = parms;
            this.model = model;
            this.representation = representation;
            this.filter = filter;
            
            if (this.parms.statistics)
                % store statistics of the representation
                this.results = struct(...
                    'f_relErrorMean', [],...
                    'f_relErrorMedian', [],...
                    'f_mean', [],...
                    'f_var', [],...
                    'f_corrcoef', [],...
                    'f_ensSpectrum', [],...
                    'f_rmse', [],...
                    'f_summary', [],...
                    'a_relErrorMean', [],...
                    'a_relErrorMedian', [],...
                    'a_mean', [],...
                    'a_var', [],...
                    'a_corrcoef', [],...
                    'a_ensSpectrum', [],...
                    'a_rmse', [],...
                    'a_summary', [],...
                    'truth', []...
                    );
                
                if (this.parms.stat_relErrors)
                    this.results.f_relErrorMean = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    this.results.a_relErrorMean = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    this.results.f_rmse = zeros(this.model.timeStepCount(), 1);
                    this.results.a_rmse = zeros(this.model.timeStepCount(), 1);
                    this.results.f_corrcoef = zeros(this.model.timeStepCount(), 1);
                    this.results.a_corrcoef = zeros(this.model.timeStepCount(), 1);
                    
                    if (this.parms.stat_summary)
                        this.results.f_relErrorMedian = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                        this.results.a_relErrorMedian = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    end
                end
                
                if (this.parms.stat_mean)
                    this.results.f_mean = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    this.results.a_mean = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                end
                
                if (this.parms.stat_var)
                    this.results.f_var = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    this.results.a_var = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                end
                
                if (this.parms.stat_summary)
                    this.results.f_summary = zeros(this.model.timeStepCount(), this.model.deterministicDimension, 5);
                    this.results.a_summary = zeros(this.model.timeStepCount(), this.model.deterministicDimension, 5);
                end
                
                if (this.parms.stat_skewness)
                    this.results.f_skewness = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    this.results.a_skewness = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                end
                
                if (this.parms.stat_kurtosis)
                    this.results.f_kurtosis = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                    this.results.a_kurtosis = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
                end
                
                if (this.parms.stat_ksdensity)
                    this.results.f_ksdensity_x = zeros(this.model.timeStepCount(), this.model.deterministicDimension, this.ksdensityValues);
                    this.results.f_ksdensity_p = zeros(this.model.timeStepCount(), this.model.deterministicDimension, this.ksdensityValues);
                    this.results.a_ksdensity_x = zeros(this.model.timeStepCount(), this.model.deterministicDimension, this.ksdensityValues);
                    this.results.a_ksdensity_p = zeros(this.model.timeStepCount(), this.model.deterministicDimension, this.ksdensityValues);
                end
                
                if (this.parms.stat_spectrum && isa(representation, 'representations.Ensemble'))
                    this.results.f_ensSpectrum = zeros(this.model.timeStepCount(), min(this.representation.sampleSize, this.model.deterministicDimension));
                    this.results.a_ensSpectrum = zeros(this.model.timeStepCount(), min(this.representation.sampleSize, this.model.deterministicDimension));
                end
                
                
            else
                % store a sample of the representation so that we can later
                % compute statistics or any functional
                this.results = struct(...
                    'f_sample', [],...
                    'a_sample', [],...
                    'truth', []...
                    );
                
                this.results.f_sample = zeros(this.model.timeStepCount(), this.model.deterministicDimension, this.representation.sampleSize);
                this.results.a_sample = zeros(this.model.timeStepCount(), this.model.deterministicDimension, this.representation.sampleSize);
            end
            
            if (this.parms.stat_truth)
                this.results.truth = zeros(this.model.timeStepCount(), this.model.deterministicDimension);
            end
            
        end
        
        function storeForecast(this, t)
            if (this.parms.statistics)
                m = this.representation.mean();
                ds = this.model.deterministicDimension;
                
                if (this.parms.stat_relErrors)
                    this.results.f_relErrorMean(t,:) = sqrt((m(1:ds) - this.model.truth).^2)./sqrt((this.model.truth).^2);
                    this.results.f_rmse(t,:) = this.representation.rmse(this.model.truth);
                    
                    if (ds > 1)
                        C = corrcoef([this.model.truth m(1:ds)]);
                        this.results.f_corrcoef(t,:) = C(1,2);
                    end
                    
                    if (this.parms.stat_summary)
                        s = this.representation.summary();
                        med = s(:,3);
                        this.results.f_relErrorMedian(t,:) = sqrt((med(1:ds) - this.model.truth).^2)./sqrt((this.model.truth).^2);
                    end
                end
                
                if (this.parms.stat_mean)
                    this.results.f_mean(t,:) = m;
                end
                
                if (this.parms.stat_var)
                    this.results.f_var(t,:) = this.representation.var();
                end
                
                if (this.parms.stat_skewness)
                    this.results.f_skewness(t,:) = this.representation.skew();
                end                
                
                if (this.parms.stat_kurtosis)
                    this.results.f_kurtosis(t,:) = this.representation.kurt();
                end
                
                if (this.parms.stat_ksdensity)
                    ens = this.representation.getEnsemble();
                    
                    for i=1:this.model.deterministicDimension
                        try
                            [~,p,x] = tools.kde(ens(i,:), this.ksdensityValues);
                            this.results.f_ksdensity_p(t,i,:) = p;
                            this.results.f_ksdensity_x(t,i,:) = x;
                        catch exception
                            if (strcmp(exception.identifier,'MATLAB:dimagree'))
                                warning('MDAPACK:Statistics:KSDensity', 'could not compute ksdensity estimate');
                                % Otherwise, just let the error propagate.
                            else
                                throw(exception); 
                            end
                        end
                    end
                end

                if (this.parms.stat_summary)
                    if ~exist('s','var')
                        s = this.representation.summary();
                    end
                    this.results.f_summary(t,:,:) = s;
                end
                
                if (this.parms.stat_spectrum && isa(this.representation, 'representations.Ensemble'))
                    E = this.representation.ensembleM;
                    E = tools.sampling.removeBias(E); % remove ensemble mean
                    [~, S, ~] = svd(E,0);
                    this.results.f_ensSpectrum(t,:) = diag(S);
                end
            else
                % store representation
                this.results.f_sample(t,:,:) = this.representation.getEnsemble();
            end
            
            if (this.parms.stat_truth)
                this.results.truth(t,:) = this.model.truth;
            end
        end
        
        function storeAnalysis(this, t)
            
            if (this.parms.statistics)
                m = this.representation.mean();
                ds = this.model.deterministicDimension;
                
                if (this.parms.stat_relErrors)
                    this.results.a_relErrorMean(t,:) = sqrt((m(1:ds) - this.model.truth).^2)./sqrt((this.model.truth).^2);
                    this.results.a_rmse(t,:) = this.representation.rmse(this.model.truth);
                    
                    if (ds > 1)
                        C = corrcoef([this.model.truth m(1:ds)]);
                        this.results.a_corrcoef(t,:) = C(1,2);
                    end
                    
                    if (this.parms.stat_summary)
                        s = this.representation.summary();
                        med = s(:,3);
                        this.results.a_relErrorMedian(t,:) = sqrt((med(1:ds) - this.model.truth).^2)./sqrt((this.model.truth).^2);
                    end
                end
                
                if (this.parms.stat_mean)
                    this.results.a_mean(t,:) = m;
                end
                
                if (this.parms.stat_var)
                    this.results.a_var(t,:) = this.representation.var();
                end
                                
                if (this.parms.stat_skewness)
                    this.results.a_skewness(t,:) = this.representation.skew();
                end
                
                if (this.parms.stat_kurtosis)
                    this.results.a_kurtosis(t,:) = this.representation.kurt();
                end
                
                if (this.parms.stat_ksdensity)
                    ens = this.representation.getEnsemble();
                    
                    for i=1:this.model.deterministicDimension
                        try
                            [~,p,x] = tools.kde(ens(i,:), this.ksdensityValues);
                            this.results.a_ksdensity_p(t,i,:) = p;
                            this.results.a_ksdensity_x(t,i,:) = x;
                        catch exception
                            if (strcmp(exception.identifier,'MATLAB:dimagree'))
                                warning('MDAPACK:Statistics:KSDensity', 'could not compute ksdensity estimate');
                                % Otherwise, just let the error propagate.
                            else
                                throw(exception); 
                            end
                        end
                    end
                end
                
                if (this.parms.stat_summary)
                    if ~exist('s','var')
                        s = this.representation.summary();
                    end
                    this.results.a_summary(t,:,:) = s;
                end
                
                if (this.parms.stat_spectrum && isa(this.representation, 'representations.Ensemble'))
                    E = this.representation.ensembleM;
                    E = tools.sampling.removeBias(E); % remove ensemble mean
                    [~, S, ~] = svd(E,0);
                    this.results.a_ensSpectrum(t,:) = diag(S);
                end
            else
                % store sample
                this.results.a_sample(t,:,:) = this.representation.getEnsemble();
            end
        end
        
        function storeFilterStatistics(this, t)
            s = this.filter.getStatistics();
            names = fieldnames(s);
            
            for j = 1:length(names)
                if isfield(this.results, ['filter_' names{j}])
                    this.results.(['filter_' names{j}]) = cat(ndims(s.(names{j})), this.results.(['filter_' names{j}]), s.(names{j}));
                else
                    %                     this.results = setfield(this.results, names{j},
                    %                     s.(names{j}));
                    this.results.(['filter_' names{j}]) = s.(names{j});
                end
            end
        end
        
        function results = getResults(this)
            results = this.results;
        end
    end
end