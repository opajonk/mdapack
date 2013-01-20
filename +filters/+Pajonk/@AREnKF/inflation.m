function X = inflation(this, X)
if strcmp(this.opts.inflation, 'adaptive')
    error('not yet implemented');
elseif strcmp(this.opts.inflation, 'fixed') && this.opts.inflation_fixed_factor > 1.0
    m = mean(X, 2);
    X = this.opts.inflation_fixed_factor .* (X - repmat(m, 1, size(X,2))) + repmat(m, 1, size(X,2));
end
end

