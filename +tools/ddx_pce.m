function [ c_new b_new ] = ddx_pce( coefficients, basis )
%DDX_PCE Summary of this function goes here
%   Detailed explanation goes here

% check dimensions
assert(size(coefficients, 2) == size(basis, 1), 'coefficients and basis have inconsistent dimension');

alpha_max = size(coefficients, 2);
dim = size(coefficients, 1);

c_new = repmat(coefficients, [1,1,dim]);
b_new = repmat(basis, [1,1,dim]);

for alpha = 1:alpha_max
    for ddw = 1:dim
        if (b_new(alpha,ddw,ddw) == 0)
            % these we could delete for efficiency
            c_new(:,alpha,ddw) = 0.0;
            b_new(alpha,:,ddw) = 0.0;
        else
            % d(H_n(w))/dw = 2*n*H_{n-1}(w)
            n = b_new(alpha,ddw,ddw);
            b_new(alpha,ddw,ddw) = b_new(alpha,ddw,ddw) - 1;
            c_new(ddw,alpha,ddw) = 2 * n * c_new(ddw,alpha,ddw);
        end
    end
end
end

