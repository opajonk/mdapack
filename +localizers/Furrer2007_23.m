function [ y ] = Furrer2007_23( p, N, x )
% FURRER2007_23 Optimal taper function for covariance function p and ensemble size N, evaluated at x
%
% FURRER2007_23(p, N, x) computes an "optimal" taper function (according to Furrer2007, section 4.3)
%   from a covariance function p. The only parameter is the ensemble size N. It is evaluated at x, where
%   x may be a vector quantity. Then it is evaluated for each vector entry.

y = 1 ./ ( 1 + ( 1 + p(0)^2 ./ p(x).^2)./N );
end

