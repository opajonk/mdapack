function [ wX l] = wtrans( this, X )
[c, l] = wavedec(X(:,1), this.opts.transformation_wavelet_level, this.Lo_D, this.Hi_D);

wX = zeros(length(c), size(X,2));
for i = 1:size(X,2)
    [c, ~] = wavedec(X(:,i), this.opts.transformation_wavelet_level, this.Lo_D, this.Hi_D);
    wX(:,i) = c;
end

assert(isempty(find(~logical(isfinite(wX)), 1)), 'MDAPACK:invalidRepresentation', 'the TRANSFORMED PRIOR ENSEMBLE is invalid as it contains NaN or Inf values - we have to abort');
end

