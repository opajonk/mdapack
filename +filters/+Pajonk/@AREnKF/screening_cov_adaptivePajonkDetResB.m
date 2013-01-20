function [L mL] = screening_cov_adaptivePajonkDetResB(this, A, HA, H, model) %#ok<INUSD>
% The difference between A and B is that B uses the expected covariance
% (the full order estimate) as mean, whereas A uses the rank-1 sample mean.

loss = @(x,y) (x.^2)./(y.^2); % quadratic loss
% loss = @(x,y) (x)./(y); % linear loss

% first compute the ensemble localiser
C = zeros(this.n,this.m,this.N-1);
C(:,:,1) = (A(:,2:end) * HA(:,2:end)') / (this.N - 2);
for i=2:this.N-2
    C(:,:,i) = (A(:,[1:i-1,i+1:end]) * HA(:,[1:i-1,i+1:end])') / (this.N - 2);
end
C(:,:,this.N-1) = (A(:,1:end-1) * HA(:,1:end-1)') / (this.N - 2);

stdC = std(C,1,3);
meanC = (A*HA') / (this.N - 1); % HERE IS THE DIFFERENCE (1/2)

r = loss(stdC,meanC);

L = (1-r./(this.N-1))./(r+1);

L(L<0) = 0.0;


% then compute the measurement localiser
mC = zeros(this.m,this.m,this.N-1);
mC(:,:,1) = (HA(:,2:end) * HA(:,2:end)') / (this.N - 2);
for i=2:this.N-2
    mC(:,:,i) = (HA(:,[1:i-1,i+1:end]) * HA(:,[1:i-1,i+1:end])') / (this.N - 2);
end
mC(:,:,this.N-1) = (HA(:,1:end-1) * HA(:,1:end-1)') / (this.N - 2);

stdmC = std(mC,1,3);
meanmC = (HA*HA') / (this.N - 1); % HERE IS THE DIFFERENCE (2/2)

r = loss(stdmC,meanmC);

mL = (1-r./(this.N-1))./(r+1);

mL(mL<0) = 0.0;
end

