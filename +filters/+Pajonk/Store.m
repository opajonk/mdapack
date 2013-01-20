classdef Store < filters.Filter
    %Store filter (simply stores the representation to files)
    
    properties
        opts;
        num;
    end
    
    methods
        function this = Store(varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addOptional('filenamePattern', 'filters.Store_%d.mat', ...
                @(x)any(strcmpi(x,{'default','advanced'})));
            
            p.parse(varargin{:});
            this.opts = p.Results;
            
            this.num = 0;
        end
        
        function this = update(this, model, representation) %#ok<INUSD>
            this.num = this.num + 1;
            save(sprintf(this.opts.filenamePattern, this.num), '-v7.3', 'representation', 'model', 'this');
        end
        
        function str = char(this)
            str = class(this);
        end
    end
end

