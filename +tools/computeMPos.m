% This method locates the maxima of the measurement operator in
% the model domain. This we use to construct a localizer for
% the measurements. Otherwise see Fertig2007a, appendix A&B.
function mpos = computeMPos(measurementDimension, deterministicDimension, h)
mresMax = ones(measurementDimension, 1)* -inf;
mpos = zeros(measurementDimension, 1);

for k = 1:deterministicDimension
    test = zeros(deterministicDimension,1);
    test(k,1) = 1;
    mres = h(test);
    
    for j = 1:measurementDimension
        if mres(j) > mresMax(j)
            mpos(j) = k;
        end
    end
    
    mresMax = max(mres, mresMax);
end
end