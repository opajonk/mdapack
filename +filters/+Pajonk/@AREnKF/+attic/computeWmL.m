function [mL] = computeWmL(this,Lo_D,Hi_D,HX,model,pos)
% obtain memory for the operator
waveletSpaceMeasurementSize = size(HX,1);

% set up some variables
T = 0.01;
M = model.measurementDimension;
h = model.measureOp();
H = h(eye(model.deterministicDimension));
H = H*H';

% decompose the measurement operator into wavelet space
wHt = zeros(waveletSpaceMeasurementSize,M);
for i = 1:M
    [wHt(:,i) ~] = wavedec(full(H(:,i)),this.waveletLevel,Lo_D,Hi_D);
end
%             % compute wavelet-space localizers if necessary
%             if strcmp(this.waveletTransform, 'measurement') || strcmp(this.waveletTransform, 'both')
%                 mL = computeWmL(this,Lo_D,Hi_D,HX,model);
%             end
%
%             if strcmp(this.waveletTransform, 'grid') || strcmp(this.waveletTransform, 'both')
%                 L = computeWL(this,Lo_D,Hi_D,X,HX,model);
%             end

% decompose the measurement operator into wavelet space
wH = zeros(waveletSpaceMeasurementSize,waveletSpaceMeasurementSize);
for i = 1:waveletSpaceMeasurementSize
    [wH(i,:) ~] = wavedec(full(wHt(i,:)),this.waveletLevel,Lo_D,Hi_D);
end

% do cut-off
wH = wH(pos:end,pos:end);
mL = conv2(wH,ones(2),'same');
mL = conv2(mL,ones(2),'same');
return;

% all coefficients above threshold mean "enough influence to be
% considered"
mL = zeros(size(wH));
mL(abs(wH)>T) = 1;

end

