function [L mL] = screening_cov_adaptiveZhang2010aBootstrap(this, A, HA, H, model) %#ok<INUSD>
loss = @(x,y) (x.^2)./(y.^2); % quadratic loss
% loss = @(x,y) (x)./(y); % linear loss

bStap = 1000;

% first compute the ensemble localiser
C = zeros(this.n,this.m,bStap);

for i=1:bStap
    C(:,:,i) = A(:,round(1+rand(1,this.N)*(this.N-1))) * HA(:,round(1+rand(1,this.N)*(this.N-1)))';
end

stdC = std(C,1,3);
meanC = A*HA';

r = loss(stdC,meanC);

L = (1-r./(this.N-1))./(r+1);

L(L<0) = 0.0;


% then compute the measurement localiser
mC = zeros(this.m,this.m,bStap);
for i=1:bStap
    mC(:,:,i) = HA(:,round(1+rand(1,this.N)*(this.N-1))) * HA(:,round(1+rand(1,this.N)*(this.N-1)))';
end

stdmC = std(mC,1,3);
meanmC = HA*HA';

r = loss(stdmC,meanmC);

mL = (1-r./(this.N-1))./(r+1);

mL(mL<0) = 0.0;
end

