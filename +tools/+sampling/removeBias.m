function [ S ] = removeBias( S )
%REMOVEBIAS Remove any column mean (bias) from a sample

% S = S - repmat(mean(S, 2), 1, size(S,2));
S = bsxfun(@minus, S, mean(S,2));
end

