function [ x ] = customMeasurementOperator(representation, H)
%CUSTOMMEASUREMENTOPERATOR Summary of this function goes here
%   Detailed explanation goes here
if (isa(representation, 'numeric')) % this is the case where we measure the "truth"
    x = H*representation;
elseif (isa(representation, 'representations.Ensemble'))
    x = H*representation.ensembleM;
elseif (isa(representation, 'representations.PCE'))
    x = H*representation.coefficients;
else
    error('customMeasurementOperator(): unsupported representation');
end

end

