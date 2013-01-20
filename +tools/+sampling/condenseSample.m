function [ sample ] = condenseSample( sourceSample, targetSize )
%CondenseSample Advanced sampling according to Evensen2004, sect. 4
assert(isfloat(sourceSample));
assert(isnumeric(targetSize));
targetSize = round(targetSize);
assert(targetSize > 1);
sourceSize = size(sourceSample, 2);
assert(targetSize < sourceSize);

% Remove any possible mean from the sample. We add it back later.
meanSourceSample = mean(sourceSample,2);
for i = 1:sourceSize
    sourceSample(:,i) = sourceSample(:,i) - meanSourceSample;
end

% h = figure(234);
% plot(std(sourceSample,0,2), 'g-');
% hold on;
% plot(mean(sourceSample,2), 'r-');

% SVD of the source sample and a random matrix which we need to create the
% smaller sample in the next step.
[U, SIGMA, ~] = svd(sourceSample, 0);
[~, ~,    V1] = svd(randn(targetSize));

beta = sourceSize / targetSize;

% Create the smaller sample using only the first targetSize singular
% vectors and singular values of sourceSample. As they are sorted by
% magnitude of singular value we capture the more important singular
% vectors with this proceduce. Scale by 1/sqrt(beta) to get correct
% variance of sourceSample back into sample.
sample = U(:,1:min(targetSize, size(U,2))) * (SIGMA(1:min(targetSize, size(SIGMA,1)),1:targetSize) .* 1/sqrt(beta)) * V1'; % scale to get correct variance

% Remove any mean that may have been introduced by this procedure and add
% back the mean of meanSourceSample.
meanSample = mean(sample,2);
for i = 1:targetSize
    sample(:,i) = sample(:,i) + meanSourceSample - meanSample;
end

% plot(std(sample,0,2), 'b-');
% plot(mean(sample,2), 'y-');
% pause;
% close(h);

end

