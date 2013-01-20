function this = update(this, model, representation)
[X sample] = representation.getEnsemble();

% check validity of ensemble
assert(isempty(find(~logical(isfinite(X)), 1)), 'MDAPACK:invalidRepresentation', 'the PRIOR ENSEMBLE is invalid as it contains NaN or Inf values - we have to abort');

% Obtain some helper variables and the measurements,
% from the model.
h = model.measureOp();
HX = h(representation);
this.N = size(X, 2);
this.n = size(X, 1);
this.m = size(HX, 1);

assert(~this.opts.normalization_fourier_enabled || this.n > 1,'MDAPACK:invalidArgument','using fourier transform with a scalar model is not clever');

% really simple inflation of the forecast ensemble according to
% Evensen2009a, Eq. (81)
if strcmp(this.opts.inflation, 'simple')
    mX = mean(X, 2);
    X = this.opts.inflation_simple_factor .* bsxfun(@plus, 1.01 * bsxfun(@minus, X, mX), mX);
end

% This method locates the maximum of the measurement operator and
% "assigns" the measurement location to that maximum in
% the model domain. In case of point measurements this is the
% same as directly using the measurement operator.
% The resulting H we use to construct a localizer for
% the measurements, mL. Another possibility could be to use the
% "support". See also Fertig2007a, appendix A&B.
H = zeros(model.measurementDimension, model.deterministicDimension);
mpos = tools.computeMPos(model.measurementDimension, model.deterministicDimension, model.measureOp());
for i = 1:model.measurementDimension
    H(i,mpos(i)) = 1;
end

if (this.opts.normalization_gridPoint_enabled)
    % Build "gridpoint normalization" which removes the smallest grid
    % point heterogeneities. Such we obtain a correlation matrix
    % later on, e.g. when computing A*HA', instead of a covariance
    % matrix. See Deckmyn2005, p. 1283.
    A = bsxfun(@minus, X, mean(X, 2));
    HA = bsxfun(@minus, HX, mean(HX, 2));
    B = spdiags(std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
    iB = spdiags(1./std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
    iHB = spdiags(1./std(HA,0,2),0,model.measurementDimension,model.measurementDimension);
    
    X = iB*X;
    HX = iHB*HX;
end

if (this.opts.normalization_fourier_enabled)
    % perform a Fourier space normalization -- this means in
    % particular that the resulting average variance spectrum
    % is white.
    fX = fft(X);
%     HXs = fft(HXs);
    A = bsxfun(@minus, fX, mean(fX, 2));
%     HA = bsxfun(@minus, HXs, mean(HXs, 2));
    F = spdiags(std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
    iF = spdiags(1./std(A,0,2),0,model.deterministicDimension,model.deterministicDimension);
%     HF = spdiags(std(HA,0,2),0,model.deterministicDimension,model.measurementDimension);
%     iHF = spdiags(1./std(HA,0,2),0,model.deterministicDimension,model.measurementDimension);
    X = ifft(full(iF*fX));
%     HXs = ifft(iHF*HXs);
end

% compute the grid point anomalies in any case so we can "have
% a look" at them and compare
%            HAs = bsxfun(@minus, HXs, mean(HXs, 2)); %#ok<NASGU>
%            As = bsxfun(@minus, Xs, mean(Xs, 2)); %#ok<NASGU>



% obtain the evidence and the evidence covariance
d = model.measure();
R = model.measureCov();

% optionally perform gridPoint normalization
if (this.opts.normalization_gridPoint_enabled)
    d = iHB*d;
    Rs = full(chol(iHB*model.measureCov()*iHB'));
    R = Rs*Rs';
end



% sample from measurement error to create measurement ensemble
if strcmp(this.opts.sampling_mode, 'advanced') && this.opts.sampling_advanced_factor > 1.0
    % Advanced sampling
    if (this.opts.normalization_gridPoint_enabled)
        E = iHB*model.measureErr(this.opts.sampling_advanced_factor*this.N, this.measurementEnsembleRng);
    else
        E = model.measureErr(this.opts.sampling_advanced_factor*this.N, this.measurementEnsembleRng);
    end
    E = tools.sampling.condenseSample(E, this.N);
else
    % Default sampling
    if (this.opts.normalization_gridPoint_enabled)
        E = iHB*model.measureErr(this.N, this.measurementEnsembleRng);
    else
        E = model.measureErr(this.N, this.measurementEnsembleRng);
    end
end

if this.opts.sampling_order >= 1
    % First order correct sampling for measurement anomalies/deviations/...
    E = tools.sampling.removeBias(E);
    
    if this.opts.sampling_order >= 2
        % Second order correct sampling for measurement anomalies/deviations/...
        E = bsxfun(@times, E, sqrt(diag(R))./sqrt(var(E,0,2)));
    end
end

% optionally perform fourier normalization on the evidence covariance
% there, we need to use the ensemble representation from Es
% if (this.opts.normalization_fourier_enabled)
%     Es = fft(Es);
%     Es = ifft(iHF * Es);
%     
%     Rs = full(chol(Es*Es'/(this.N-1)));
% end

% create measurement ensemble
D = bsxfun(@plus, E, d);


% sanity check, performed in any case to allow making on the fly wavelet
% transforms
if (wmaxlev(this.n, this.opts.transformation_wavelet_type) < this.opts.transformation_wavelet_level)
    this.opts.transformation_wavelet_level = wmaxlev(this.n, this.opts.transformation_wavelet_type);
end
    

% optionally perform the wavelet transform
if (this.opts.transformation_wavelet_enabled)
    if (any(strcmp(this.opts.transformation_wavelet_mode, {'grid','both'})))
        [X this.wlX] = this.wtrans(X);
    end
    
    if (any(strcmp(this.opts.transformation_wavelet_mode, {'measurement','both'})))
        [HX this.wlHX] = this.wtrans(HX);
    end
    
    if(this.opts.transformation_wavelet_localization_enabled)
        this.prepare_wL_MinMax(model);
    end
end

% compute ensemble deviations/anomalies/whatever-they-call-this...
A = bsxfun(@minus, X, mean(X, 2));
HA = bsxfun(@minus, HX, mean(HX,2));

% optionally compute localizers
if (~strcmp(this.opts.localization, 'none') && isempty(this.covLoc))
    this.prepare_localization(model,H);
end

% optionally compute screeners
if (strcmp(this.opts.screening, 'cov'))
    [S mS] = feval(this.opts.screening_function_cov, this, A, HA, H, model);
end

P = (HA * HA');

assert(isempty(find(~logical(isfinite(P)), 1)), 'MDAPACK:invalidRepresentation', 'the VAR(HX) ESTIMATE is invalid as it contains NaN or Inf values - we have to abort');

% optionally perform localization of measurement error covariance
if (any(strcmp(this.opts.localization, {'measurement','both'})))
    P = this.mCovLoc .* P;
end

% optionally perform screening of measurement error covariance
if (strcmp(this.opts.screening, 'cov'))
    P = mS .* P;
end

P = (this.N-1) * R + P;

C = A * HA';

assert(isempty(find(~logical(isfinite(C)), 1)), 'MDAPACK:invalidRepresentation', 'the COV(C,HX) ESTIMATE is invalid as it contains NaN or Inf values - we have to abort');

% optionally perform localization of measurement-to-grid error covariance
if (any(strcmp(this.opts.localization, {'measurement','both'})))
    C = this.covLoc .* C;
end

if (strcmp(this.opts.screening, 'cov'))
    C = S .* C;
end


%% use one of the different solution methods to compute Kalman gain and update
if strcmp(this.opts.solving, 'MATLAB')
    
    X = this.inflation(X);
    
    % use the simple MATLAB inversion
    X = X + C * (P \ (D - HX));
    
elseif strcmp(this.opts.solving, 'eigen')
    % solve by truncated eigen decomposition
    
    [eigv, eigd] = eig(P);
    
    n_eig = length(d);
%     n_eig = min(this.N,length(d));
    
    weights = eigv(1:n_eig,1:n_eig) * diag(1./diag(eigd(1:n_eig,1:n_eig))) * eigv(1:n_eig,1:n_eig)';
    K = C * weights;
    
    if (this.opts.transformation_wavelet_localization_enabled)
        K = this.waveletGainLoc .* K;
    end
    
    if (strcmp(this.opts.screening, 'gain'))
        K = feval(this.opts.screening_function_gain, this, K, A, HA, R, model);
    end
    
    X = this.inflation(X);
    
    X = X + K * (D - HX);
    
elseif strcmp(this.opts.solving, 'svd')
    % This is the "alternative solution for large measurement dimensions"
    % like [Evensen2003]. It also has the advantage of not inverting the
    % normal equation HA*HA', as this suares the condition number of that
    % matrix and unnecessarily amplifies noise.
    
    % compute (HA*HA' + E*E')^-1 = ((HA + E)(HA + E)')^-1 =
    % [compute SVD]
    % = ((U*S*V')(U*S*V')')^-1 = (U*S*V'*V*S'*U')^-1 = (U*S*V'*V*S'*U')^-1
    % = (U*S*S'*U')^-1 = U*(S*S')^-1*U'
    [U S ~] = svd(HA + E,'econ');
    
    s = diag(S);
    invss = 1./(s.^2);
    
    % compute how many singular values are retained in the pseudo inversion
    tol = max(this.m,this.N) * eps(max(s));
    r = sum(s > tol); % only the ones which are bigger than a certain threshold
    
    weights = U(:,1:r) * diag(invss(1:r)) * U(:,1:r)';
    
    K = C * weights;
    
    if (this.opts.transformation_wavelet_localization_enabled)
        K = this.waveletGainLoc .* K;
    end
    
    if (strcmp(this.opts.screening, 'gain'))
        K = feval(this.opts.screening_function_gain, this, K, A, HA, R, model);
    end
    
    X = this.inflation(X);
    
    X = X + K * (D - HX);
else
    error('unknown solution method - this should never happen and is a programming error');
end

assert(isempty(find(~logical(isfinite(X)), 1)), 'MDAPACK:invalidRepresentation', 'the POSTERIOR ENSEMBLE is invalid as it contains NaN or Inf values - we have to abort');

% optionally transform back from wavelet space
if (this.opts.transformation_wavelet_enabled) && (any(strcmp(this.opts.transformation_wavelet_mode, {'grid','both'})))
    X = this.iwtrans(X,this.wlX);

    assert(isempty(find(~logical(isfinite(X)), 1)), 'MDAPACK:invalidRepresentation', 'the TRANSFORMED POSTERIOR ENSEMBLE is invalid as it contains NaN or Inf values - we have to abort');
end

% optionally transform back from normalized fourier space
if (this.opts.normalization_fourier_enabled)
    fX = fft(X);
    X = ifft(F*fX);
end

% optionally transform back into grid point space if necessary
if (this.opts.normalization_gridPoint_enabled)
    X = B*X;
end

representation.setEnsemble(X, sample);
end

