function [L mL] = screening_cov_adaptivePajonkDetResC(this, A, HA, H, model) %#ok<INUSD>
% The difference between A and B is that B uses the expected covariance
% (the full order estimate) as mean, whereas A uses the rank N-1 sample mean.

loss = @(x,y) (x.^2)./(y.^2); % quadratic loss
% loss = @(x,y) (x)./(y); % linear loss

sampleSize = this.N;
maxBootstrapSize = 10;

varC = zeros(this.n,this.m);
meanC = A(:,round(1+rand(1,sampleSize)*(this.N-1))) * HA(:,round(1+rand(1,sampleSize)*(this.N-1)))' ./ (sampleSize-1);

while true
    %     disp(num2str(norm(meanCnew - meanC)/norm(meanCnew)));
    
    newC = A(:,round(1+rand(1,sampleSize)*(this.N-1))) * HA(:,round(1+rand(1,sampleSize)*(this.N-1)))' ./ (sampleSize-1);
    
    meanCnew = meanC + (newC - meanC)/(counter+1);
    varC = (1 - 1/counter) .* varC + (counter+1).*(meanCnew - meanC).^2;
    
    counter = counter + 1;
    
    if norm(meanCnew - meanC)/norm(meanC) < 0.01 || counter > maxBootstrapSize
        break;
    end
    
    meanC = meanCnew;
end
meanC = meanCnew;

r = (varC)./(meanC.^2);

% L = (sampleSize-r)./( (sampleSize-1).*r + sampleSize);
L = (1-r./(sampleSize-1))./(r+1);

L(L<0) = 0.0;

mL = H*L;
end

