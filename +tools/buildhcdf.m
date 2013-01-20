function [ hcdf nbins ] = buildhcdf( x )
N = size(x,2);
nbins = floor(sqrt(N));


if size(x,1) > 1
    hcdf = struct([]);
    
    idx_lower=1;
    idx_lower_previous=1;
    j=1;
    xs = sortrows(x',1)';
    
    while (idx_lower<N)
        idx_upper =min(idx_lower + floor(nbins/3*2),N); % skip at least half the size of a bin, or to the end
        while(idx_upper < N && xs(1,idx_upper) == xs(1,(idx_upper)+1)) % search for the next different value
            idx_upper = idx_upper+1;
        end
        % a single-element last class is not good
        if (idx_upper+1 == N)
            idx_upper = idx_upper+1;
        end
        
        if (min(xs(1,idx_lower:idx_upper)) == max(xs(1,idx_lower:idx_upper)))
            % class contains only the same value - not so good
            if (idx_upper < N)
                % we can add one more value *phew*
                idx_upper = idx_upper + 1;
            else
                % kill this class and join with previous ones
                if j > 1
                    hcdf(j-1).sup = xs(1,N);
                    hcdf(j-1).cdf = tools.buildhcdf(xs(2:end,idx_lower_previous:N));
                    idx_lower = N;
                    continue;
                else
                    % ok, then we need to stay with this one-valued
                    % class...
                end
            end
        end
        hcdf(j).min = xs(1,idx_lower);
        hcdf(j).sup = xs(1,min(idx_upper+1,N));
        hcdf(j).cdf = tools.buildhcdf(xs(2:end,idx_lower:idx_upper));
        
        idx_lower_previous = idx_lower;
        idx_lower = idx_upper+1;
        
        j = j+1;
    end
    
    nbins = length(hcdf);
else
%     if min(x)==max(x)
%         hcdf = @(c) c<=x(1);
%     else
N = length(x);

if (N > 20)
    % compute piecewise hermite interpolation function
    [~,~,xmesh,cdf] = tools.kde(x);
    cdf_poly = pchip(xmesh, cdf);
    hcdf = @(c) ppval(c, cdf_poly);
    nbins = 1;
else
    
    % use the erf method
%     x = sort(x);
%     intv = max(x) - min(x);
    hcdf = @(c) arrayfun(@(d) mean(arrayfun(@(r) randstep(d,r),x)), c);
%     hcdf = @(c) arrayfun(@(d) mean(erf((d - x)*(intv*1000000)) + 1)*0.5, c);
end
%     end
end
end


function v = randstep(d,x)
if randi(2,1) == 1, v=d>x; else v=d>=x; end
end