function [] = prepare_wL_MinMax(this, model)

if (isempty(this.waveletGainLoc))
    discretisationSteps = 100;
    dist = model.distanceMatrix();
    
    sigma_max = model.decorrelationLength();
    sigma_min = model.measurementResolution();
    
    this.waveletGainLoc = zeros(this.n, this.m);
    
    % determine relative importance
    for i = sigma_min:(sigma_max-sigma_min)/discretisationSteps:sigma_max
        discreteLocFun = this.opts.localization_function(dist, i, this.N);
        relativeImportance = abs(this.wtrans(discreteLocFun));
        this.waveletGainLoc = max(this.waveletGainLoc, relativeImportance);
    end
    
    % normalize each level to max=1.0
    for k = 1:this.m
        runningSum = 1;
        for j = 1:length(this.wlX)-1
            this.waveletGainLoc(runningSum:runningSum+this.wlX(j)-1,k) = this.waveletGainLoc(runningSum:runningSum+this.wlX(j)-1,k)./max(this.waveletGainLoc(runningSum:runningSum+this.wlX(j)-1,k));
            runningSum = runningSum + this.wlX(j);
        end
    end
    
    % perform thresholding if necessary
    if (this.opts.transformation_wavelet_localization_threshold_enabled)
        this.waveletGainLoc(this.waveletGainLoc >= this.opts.transformation_wavelet_localization_threshold) = 1;
        this.waveletGainLoc(this.waveletGainLoc < this.opts.transformation_wavelet_localization_threshold) = 0;
    end
end
end

