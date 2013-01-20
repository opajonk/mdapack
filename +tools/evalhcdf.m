function [ y ] = evalhcdf( hcdf, x)
nbins = length(hcdf);

if (nbins == 1 && isstruct(hcdf))
    y = tools.evalhcdf( hcdf(1).cdf, x(2:end));
    return;
end

if (isstruct(hcdf))
    % seach correct interval and call recursively
    
    for k = 1:nbins
        if (k==1 && x(1) < hcdf(k).sup) ... % less than minimum in first interval is OK
                || (k==nbins && x(1) >= hcdf(k).min) ... % more than maximum in last interval is OK
                || (x(1) >= hcdf(k).min && x(1) < hcdf(k).sup) % default correct interval
            y = tools.evalhcdf( hcdf(k).cdf, x(2:end));
            break;
        end
    end
    
    
%     k = a + floor(nbins/2)-1;
%     if (x(1) < hcdf(k).min);
%         y = evalhcdf(hcdf, x, a, k-1 );
%     elseif (x(1) >= hcdf(k).sup && nbins > 1);
%         y = evalhcdf( hcdf, x, min(k+1,b), b );
%     else
%         % were in the right interval
%         
%     end
    
elseif isa(hcdf, 'function_handle')
    [y] = hcdf(x);
%     if (y >= 1.0) 
%         y = 1-eps;
%     end
end
end

