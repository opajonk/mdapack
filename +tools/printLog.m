function printLog( fileName, formatString, varargin )
%PRINTLOG Print a log message
fprintf(1, '%12s: ', fileName);
fprintf(1, formatString, varargin{:});
end

