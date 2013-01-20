function [ x ] = quadraticMeasurementOperator(representation, dims)

if (isa(representation, 'numeric')) % this is the case where we measure the "truth"
    pos = 1;
    x = zeros(length(dims),1);
    for i = dims
        x(pos,:) = representation(i,:)*representation(i,:);
        pos = pos + 1;
    end
elseif (isa(representation, 'representations.Ensemble'))
    pos = 1;
    x = zeros(length(dims), size(representation.ensembleM, 2));
    for i = dims
        x(pos,:) = representation.ensembleM(i,:).*representation.ensembleM(i,:);
        pos = pos + 1;
    end
elseif (isa(representation, 'representations.PCE'))
    pos = 1;
    x = zeros(length(dims), size(representation.basis, 1));
    for i = dims
        x(pos,:) = pce_multiply(representation.coefficients(i,:), representation.basis, representation.coefficients(i,:), representation.basis, representation.basis, 'M', representation.multTensor);
        pos = pos + 1;
    end
else
    error('customMeasurementOperator(): unsupported representation');
end

end

