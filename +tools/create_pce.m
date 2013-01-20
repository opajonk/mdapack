function pce = create_pce(dim, order, vals, samples, weights)

assert(sum(size(samples) == size(vals))==2, 'vals and samples has to have same dimensions');

N = size(vals, 2);

% default weights: all equal
if nargin < 5
    weights = ones(1, N) ./ N;
end

indx_glob=multiindex(dim, order);
herm_poly = hermite_mat(order);
multivar_normalize = multiindex_factorial(indx_glob);
multivar_order = size(indx_glob,1);

k_coef = zeros(dim,length(multivar_normalize),N);
for i = 1:N
    %  compute the values for each hermite polynom for the current sample
    herm_vals = zeros(order+1, dim);
    
    for p=1:order+1
        herm_vals(p,:) = polyval(herm_poly(p,:), samples(:,i)');
    end
    
    % copy to positions
    multi_vals = zeros(length(multivar_normalize), dim);
    
    for k = 1:multivar_order
        for j = 1:dim
            multi_vals(k,j) = herm_vals(indx_glob(k,j)+1, j);
        end
    end
    
    % compute the normalized multivariate values
    psi_vals = prod(multi_vals, 2) ./ multivar_normalize;
    
    k_coef(:,:,i) = reshape(kron(psi_vals, vals(:,i)), dim, []);
    
    k_coef(:,:,i) = weights(i)*k_coef(:,:,i);
end

pce = sum(k_coef,3);