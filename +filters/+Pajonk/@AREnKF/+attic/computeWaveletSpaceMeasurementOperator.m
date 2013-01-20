function [ HW ] = computeWaveletSpaceMeasurementOperator( his,Lo_D,Hi_D,H )
%  This function can be used to optionally transform the linear
% measurement operator H to Wavelet space.
H = full(H)';
[coefficients, ~] = wavedec(H(:,1),this.waveletLevel,Lo_D,Hi_D);
HW = zeros(length(coefficients), size(H,2));
for i = 1:size(H,2)
    [coefficients, ~] = wavedec(H(:,i),this.waveletLevel,Lo_D,Hi_D);
    HW(:,i) = coefficients;
end

end

