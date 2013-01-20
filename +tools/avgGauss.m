function [ avgGaussH ] = avgGauss(nDim, mDim, sigma )
% averagingWidth = 2;
% mDim = 20;
% nDim = 40;
x = -floor(nDim/2)+1:floor(nDim/2);
avgGaussH = zeros(mDim,nDim);
for i = 1:mDim
    foo = exp(-1.0*(x).^2/(sqrt(2)*sigma).^2);
    avgGaussH(i, :) = circshift(foo,[0 floor(nDim/mDim*i + nDim/2 - 1)]);
end

end

