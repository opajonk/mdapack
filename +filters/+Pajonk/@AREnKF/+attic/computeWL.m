function [L] = computeWL(this,Lo_D,Hi_D,X,HX,model,pos)
% obtain memory for the operator
waveletSpaceDomainSize = size(X,1);
waveletSpaceMeasurementSize = size(HX,1);

% set up some variables
T = 0.01;
M = model.measurementDimension;
h = model.measureOp();
H = h(eye(model.deterministicDimension));

% decompose the measurement operator into wavelet space
wHt = zeros(waveletSpaceDomainSize,M);
for i = 1:M
    [wHt(:,i) ~] = wavedec(full(H(i,:)),this.waveletLevel,Lo_D,Hi_D);
end


if strcmp(this.waveletTransform, 'both')
    % decompose the measurement operator into wavelet space
    wH = zeros(waveletSpaceDomainSize,waveletSpaceMeasurementSize);
    for i = 1:waveletSpaceDomainSize
        [wH(i,:) ~] = wavedec(full(wHt(i,:)),this.waveletLevel,Lo_D,Hi_D);
    end
else
    wH = wHt;
end

% all coefficients above threshold mean "enough influence to be
% considered"
L = zeros(size(wH));
L(abs(wH)>T) = 1;
end

