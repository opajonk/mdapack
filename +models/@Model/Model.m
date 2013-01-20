classdef Model < handle
    %MODEL Abstract base class for all models
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        truth;
    end
    
    methods (Abstract)
        step(this, t)
        [] = createInitialEnsemble(this, ens, varargin)
        [] = createInitialPCE(this, pce, varargin)
        [err] = stepErr(this, t)
        [hasM] = hasMeasurement(this)
        [m] = measure(this, R)
        [H] = measureOp(this)
        [this] = plot(this, representation)
        [merr] = measureErr(this, N)
        [cov] = measureCov(this)
        [L] = distanceMatrix(this)
        [tCount] = timeStepCount(this)
    end
end

