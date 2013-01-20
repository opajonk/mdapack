function [this] = plot(this, representation, varargin)
%     subplot(2,1,1);
    if isempty(this.truthData)
        this.truthData = this.truth;
        this.timeScale = this.time;
    else
        this.truthData = cat(2, this.truthData, this.truth);
        this.timeScale = cat(2, this.timeScale, this.time);
    end
    
    plot(this.timeScale, this.truthData, varargin{:});
    mup = 123.0+5;
    mdown = 123.0-5;
    ylabel('random process x')
    
    axis([0, ceil(max(this.timeScale))+1, mdown, mup]);
    axis manual;
    xlabel('time');
    
    hold on;
    
    xmean = representation.mean();
    xsummary = representation.summary();
    xstd = representation.std();
    
%     step = 20/(sqrt(representation.ens_size));
%     edges = -10:step:10;
%     
%     [n,bin] = histc(representation.ensembleM(1,:), edges);
%     m = mode(bin);
%     xmean(1) = mean(edges([m, m+1]));
%     
%     [n,bin] = histc(representation.ensembleM(2,:), edges);
%     m = mode(bin);
%     xmean(2) = mean(edges([m, m+1]));
%     
%     [n,bin] = histc(representation.ensembleM(3,:), edges);
%     m = mode(bin);
%     xmean(3) = mean(edges([m, m+1]));
    
    if isempty(this.ensMeanData)
        this.ensMeanData = xmean;
    else
        this.ensMeanData = cat(2, this.ensMeanData, xmean);
    end
    if isempty(this.ensStdData)
        this.ensStdData = xstd;
    else
        this.ensStdData = cat(2, this.ensStdData, xstd);
    end
    
    if isempty(this.summaryData)
        this.summaryData = xsummary;
    else
        this.summaryData = cat(3, this.summaryData, xsummary);
    end
    
    if (size(this.summaryData, 3) == 1)
        plot(this.timeScale, this.summaryData(3), 'r-', 'LineWidth', 2);
    else
        plot(this.timeScale, squeeze(this.summaryData(1,3,:)), 'r-', 'LineWidth', 2);
    end
    %                 plot(this.timeScale, this.ensMeanData(this.plot2dDim,:), 'r-', 'LineWidth', 2);
    
    % plot quantiles p=2.5, p=25, p=50, p=75, p=97.5
    if (size(this.summaryData, 3) == 1)
        plot(this.timeScale, this.summaryData(1,1), 'r.');
        plot(this.timeScale, this.summaryData(1,2), 'r--');
        plot(this.timeScale, this.summaryData(1,4), 'r--');
        plot(this.timeScale, this.summaryData(1,5), 'r.');
    else
        plot(this.timeScale, squeeze(this.summaryData(1,1,:)), 'r.');
        plot(this.timeScale, squeeze(this.summaryData(1,2,:)), 'r--');
        plot(this.timeScale, squeeze(this.summaryData(1,4,:)), 'r--');
        plot(this.timeScale, squeeze(this.summaryData(1,5,:)), 'r.');
    end
    
    
    %                 plot(this.timeScale, this.ensMeanData(this.plot2dDim,:) + this.ensStdData(this.plot2dDim,:), 'r--');
    %                 plot(this.timeScale, this.ensMeanData(this.plot2dDim,:) - this.ensStdData(this.plot2dDim,:), 'r--');
    % plot max/min -- shows if the truth is "covered" by the ensemble
    % at all
    % plot(max(x,[],2), 'r--');
    % plot(min(x,[],2), 'r--');
    
    
    
    if(this.hasMeasurement())
        if isempty(this.mPosData)
            this.mPosData = [length(this.truthData(1,:)); this.truthData(:,end)];
%             this.mPosData = [length(this.truthData(1,:)); this.lastMeasurement];
        else
            this.mPosData = cat(2,this.mPosData, [length(this.truthData(1,:)); this.truthData(:,end)]);
%             this.mPosData = cat(2,this.mPosData, [length(this.truthData(1,:)); this.lastMeasurement]);
        end
    end
    if ~isempty(this.mPosData)
        plot(this.mPosData(1,:) * this.tStep, this.mPosData(2,:), 'ro','markerfacecolor',[0 0 0]);
    end
    
    hold off;
    
%     subplot(2,1,2);
    
    return;
    
    if isa(representation, 'representations.Ensemble')
%        hist(this.ensemble(this.plot2dDim,:));
%        [f, xi] = ksdensity(this.ensemble(this.plot2dDim,:));
%        plot(xi,f);
%        axis([mdown, mup, 0, size(this.ensemble, 2)/2 ]);
       [bw,p,x] = tools.kde( representation.ensembleM(:));
       plot(x, p, '-');
       xlim([-2 4]);
       set(gca,'YTickLabel','');
       set(gca,'YTick',[]);
%        axis([mdown, mup, 0, 2 ]);
%         empirical_density( representation.ensembleM(this.plot2dDim,:), 30, 10, 'plot_args', {'-*'} );
    elseif isa(representation, 'representations.PCE')
        sample = randn(this.deterministicDimension, 1000);
        sample = tools.sampling.removeBias(sample);
        sample = bsxfun(@times, sample, 1.0 ./sqrt(var(sample,0,2)));
        
        vals = pce_evaluate(representation.coefficients, representation.basis, sample);
        [bw,p,x] = tools.kde(vals(:));
%         p = wden(p,'heursure','s','one',2,'sym8');
%         p = p./sum(p); % normalize, so it is a PDF
        plot(x, p, '-');
        ylim([0 max(p)*1.1]);
        xlim([-3 4]);
       set(gca,'YTickLabel','');
       set(gca,'YTick',[]);
    end

end