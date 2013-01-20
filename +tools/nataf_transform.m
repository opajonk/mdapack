function [ tSamples ] = nataf_transform( samples )
%NATAF_TRANSFORM Summary of this function goes here
%   Detailed explanation goes here

tSamples = zeros(size(samples));
[N M] = size(samples);

for i=1:M
    [bandwidth,density,xmesh,cdf] = tools.kde(samples(:,i));
    
    % HACK(?): Make sure the values of the density estimate are
    % valid - otherwise xnew may contain NaNs
    % Is there a better solution for this?
    cdf(cdf<=0.0) = eps;
    cdf(cdf>=1.0) = 1-eps;
    
    % interpolate using piecewise cubic Hermite polynomials
    % these preserve monotonicity (in the way they are setup in MATLAB)
    unifx = pchip(xmesh,cdf,samples(:,i));
    tSamples(:,i) = norminv(unifx,0,1);
end


% 4) Decorrelate xnew using a square root of the inverse of the
%    correlation matrix estimated from xnew
xtilda = bsxfun(@times, bsxfun(@minus, tSamples, mean(tSamples, 1)), 1./std(tSamples));
R0 = 1/(N-1) * (xtilda' * xtilda);
iR0 = pinv(R0);
L = chol(iR0);
tSamples = L*tSamples';

end

