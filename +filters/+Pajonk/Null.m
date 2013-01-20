classdef Null < filters.Filter
    %Null filter (does nothing, e.g. for model integration runs)
    
    properties
        opts;
    end
    
    methods
        function this = Null(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            this.opts = p.Results;
        end
        
        function this = update(this, ~, ~)
            % do nothing
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

