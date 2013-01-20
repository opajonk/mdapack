function [y] = Gaspari1999_4_10(x,c)
% GASPARI1999_4_10 Compactly supported 5th order piecewise rational correlation function
%
% GASPARI1999_4_10(X,c) computes a compactly supported 5th order piecewise rational
%   correlation function. It is 1.0 at x=0.0 and smoothly drops to 0
%   for |x| > 2*c, where x is an array of distance and c is a parameter describing the support
%   width of the function. It is similar to a Gaussian and represents a homogeneous
%   and isotropic correlation function.
%
% See Gaspari1999, p. 748, "(c) Compactly supported 5th order piecewise rational functions", equation (4.10)

y = arrayfun(@(xi) scalarGaspari1999_4_10(xi,c), x);
end

function y = scalarGaspari1999_4_10(x,c)
if (0 <= abs(x) && abs(x) <= c)
    y = -1/4*(abs(x)/c)^5 + 1/2*(x/c)^4 + 5/8*(abs(x)/c)^3 - 5/3*(x/c)^2 + 1;
elseif (c < abs(x) && abs(x) <= 2*c)
    y = 1/12 * (abs(x)/c)^5 - 1/2*(x/c)^4 + 5/8*(abs(x)/c)^3 + 5/3*(x/c)^2 - 5*(abs(x)/c) + 4 - 2/3*c/abs(x);
else % x > 2*c
    y = 0;
end
end