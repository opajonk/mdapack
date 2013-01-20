% !!IMPORTANT NOTICE!!
% THIS FILE COMES FROM THE OSPREY TOOLKIT
% !!IMPORTANT NOTICE!!
%
% (including the whole toolkit was not possible due to incompatibilities;
% it has been changed in a minor way to allow for incorporation)
%
% See http://www.mathworks.com/matlabcentral/fileexchange/70-osprey

function y = percentile(x, pct)
% y = PERCENTILE(x, pct)
% Example: percentile(x, 0.10) is the value that is higher than 10%
% of the elements of x, and less than the remaining 90%.
% If the length of x is such that there is no element of x exactly
% corresponding to the 'pct' position, a weighted average of the two
% adjacent values is used.  pct must be between 0 and 1 inclusive.
%
% percentile(x, 1)   is a slow way to get max(x).  
% percentile(x, 0.5) is the median of x.
% percentile(x, 0)   is a slow way to get min(x).  
%
% If x is a matrix, percentile operates on columns, returning multiple columns.
% If pct is a vector, multiple rows are returned, one per element of pct.
%
% See also median, sort, max, min.

if (tools.nRows(x) == 1), x = x.'; end
x = sort(x);
n = ((tools.nRows(x) - 1) * pct(:) + 1);
r = rem(n, 1);
r1 = r * ones(1, tools.nCols(x));
y = (1-r1) .* x(n-r,:);
ix = find(r);			% when n=nRows(x), x(n+1,:) doesn't exist
y(ix,:) = y(ix,:) + r1(ix,:) .* x(n(ix)-r(ix)+1,:);
