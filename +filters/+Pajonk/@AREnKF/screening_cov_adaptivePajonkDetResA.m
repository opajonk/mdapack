function [L mL] = screening_cov_adaptivePajonkDetResA(this, A, HA, H, model) %#ok<INUSD>
% The difference between A and B is that B uses the expected covariance
% (the full order estimate) as mean, whereas A uses the rank-1 sample mean.

loss = @(x,y) (x.^2)./(y.^2); % quadratic loss
% loss = @(x,y) (x)./(y); % linear loss

counter = 1;

meanC = (A(:,2:end) * HA(:,2:end)') / (this.N - 2);
varC = zeros(this.n,this.m);

counter = counter + 1;

for i=2:this.N-2
    newC = (A(:,[1:i-1,i+1:end]) * HA(:,[1:i-1,i+1:end])') / (this.N - 2);
    
    meanCnew = meanC + (newC - meanC)/(counter+1);
    varC = (1 - 1/counter) .* varC + (counter+1).*(meanCnew - meanC).^2;
    
    meanC = meanCnew;
    counter = counter + 1;
end

newC = (A(:,1:end-1) * HA(:,1:end-1)') / (this.N - 2);
meanCnew = meanC + (newC - meanC)/(counter+1);
varC = (1 - 1/counter) .* varC + (counter+1).*(meanCnew - meanC).^2;
meanC = meanCnew;


r = (varC)./(meanC.^2);

L = (1-r./(this.N-1))./(r+1);

L(L<0) = 0.0;

mL = H*L;

end

