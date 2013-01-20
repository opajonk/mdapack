classdef Filter < handle
    %Filter Base class for all filters
    
    properties
    end
    
    methods
        function stat = getStatistics(this)
            stat = struct([]);
        end
    end
    
    methods (Abstract)
        this = update(this, model, representation)
        str = char(this)
    end
end