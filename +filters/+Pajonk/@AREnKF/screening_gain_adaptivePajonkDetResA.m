function [K] = screening_gain_adaptivePajonkDetResA(this, K, A, HA, R, model, H) %#ok<INUSD>

counter = 1;

meanK = (A(:,2:end) * HA(:,2:end)') / (this.N - 2) * pinv((this.N-2) * R + ( HA(:,2:end) *  HA(:,2:end)')/ (this.N - 2));
varK = zeros(size(meanK));

counter = counter + 1;

for i=2:this.N-2
    newK = (A(:,[1:i-1,i+1:end]) * HA(:,[1:i-1,i+1:end])') / (this.N - 2) * pinv((this.N-1) * R + ( HA(:,[1:i-1,i+1:end]) *  HA(:,[1:i-1,i+1:end])')/ (this.N - 2));
    
    meanKnew = meanK + (newK - meanK)/(counter+1);
    varK = (1 - 1/counter) .* varK + (counter+1).*(meanKnew - meanK).^2;
    
    meanK = meanKnew;
    counter = counter + 1;
end

newK = (A(:,1:end-1) * HA(:,1:end-1)') / (this.N - 2) * pinv((this.N-2) * R + ( HA(:,1:end-1) *  HA(:,1:end-1)')/ (this.N - 2));

meanKnew = meanK + (newK - meanK)/(counter+1);
varK = (1 - 1/counter) .* varK + (counter+1).*(meanKnew - meanK).^2;
meanK = meanKnew;


r = (varK)./(meanK.^2);

% Variant 1 of Zhang2010a
% SK = (ones(size(r))-r./(this.N-2))./(r+ones(size(r)));  
% SK(SK<0) = 0.0;

% Variant 2 of Zhang2010a, 0.6 is a variable
sigma_a = 0.035;
SK = ones(size(r))./(ones(size(r)) + r .* (ones(size(r)) ./ (sigma_a)^2)); 

% My own simple idea, with variable values
% a = 0.5;
% b = 0.5;
% SK=tanh(sqrt(r)) .* a + b; 

% K = SK .* K;

H = abs(this.wtrans_X(H'));
% H(H>0.001) = 1;


% decorrelation part
sigma = 10.0;
mDim = 20;
nDim = 40;
avgGaussH = tools.avgGauss(nDim, mDim, sigma);
dcor = abs(this.wtrans_X(avgGaussH'));
% dcor(dcor>0.001) = 1;

loc = max(dcor,H);

% loc(loc > 0.01) = 1;

loc(loc <0.01) = 0;


for k = 1:this.m
    runningSum = 1;
    for i = 1:length(this.wlX)-1
        loc(runningSum:runningSum+this.wlX(i)-1,k) = loc(runningSum:runningSum+this.wlX(i)-1,k)./max(loc(runningSum:runningSum+this.wlX(i)-1,k));
        runningSum = runningSum + this.wlX(i);
    end
end
% loc = dcor;
% loc = H;
Knonloc = K;
K = loc.* K;
end

