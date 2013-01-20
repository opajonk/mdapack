function [this] = plot(this, representation, varargin)

if(this.plotPDF)
    subplot(2,1,1);
end

if isempty(this.truthData)
    this.truthData = this.truth;
    this.timeScale = this.time;
else
    this.truthData = cat(2, this.truthData, this.truth);
    this.timeScale = cat(2, this.timeScale, this.time);
end

if strcmp(this.plotMode, '3d')
    plot3(this.truthData(1,:), this.truthData(2,:), this.truthData(3,:), varargin{:});
    axis([-20 20 -40 20 10 50]);
    axis manual;
    mup = 20;
    mdown = -20;
else
    plot(this.timeScale, this.truthData(this.plot2dDim,:), varargin{:});
    switch (this.plot2dDim)
        case 1
            mup = 20;
            mdown = -20;
            ylabel('intensity of convective motion');
        case 2
            mup = 20;
            mdown = -40;
            ylabel('temperature difference of currents');
        case 3
            mup = 50;
            mdown = 0;
            ylabel('distortion of vertical temp. profile');
        case 4
            mup = 35;
            mdown = 20;
            ylabel('parameter rho');
        case 5
            mup = 15;
            mdown = 5;
            ylabel('parameter sigma');
        case 6
            mup = 3*8/3;
            mdown = 0.3*8/3;
            ylabel('parameter beta');
    end
    
    axis([0, ceil(max(this.timeScale))+1, mdown, mup]);
    axis manual;
    xlabel('time (days)');
end

hold on;

xmean = representation.mean();
xstd = representation.std();
xsummary = representation.summary();

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

if strcmp(this.plotMode, '3d')
    % TODO How to plot the standard deviation of the ensemble in 3D?
    plot3(this.ensMeanData(1,:),this.ensMeanData(2,:),this.ensMeanData(3,:), 'r-');
else
    if (size(this.summaryData, 3) == 1)
        plot(this.timeScale, this.summaryData(this.plot2dDim,3), 'r-', 'LineWidth', 2);
    else
        plot(this.timeScale, squeeze(this.summaryData(this.plot2dDim,3,:)), 'r-', 'LineWidth', 2);
    end
    
    if (size(this.summaryData, 3) == 1)
        plot(this.timeScale, this.summaryData(this.plot2dDim,1), 'r.');
        plot(this.timeScale, this.summaryData(this.plot2dDim,2), 'r--');
        plot(this.timeScale, this.summaryData(this.plot2dDim,4), 'r--');
        plot(this.timeScale, this.summaryData(this.plot2dDim,5), 'r.');
    else
        plot(this.timeScale, squeeze(this.summaryData(this.plot2dDim,1,:)), 'r.');
        plot(this.timeScale, squeeze(this.summaryData(this.plot2dDim,2,:)), 'r--');
        plot(this.timeScale, squeeze(this.summaryData(this.plot2dDim,4,:)), 'r--');
        plot(this.timeScale, squeeze(this.summaryData(this.plot2dDim,5,:)), 'r.');
    end
end



if(this.hasMeasurement())
    if isempty(this.mPosData)
        this.mPosData = [length(this.truthData(1,:)); this.truthData(:,end)];
    else
        this.mPosData = cat(2,this.mPosData, [length(this.truthData(1,:)); this.truthData(:,end)]);
    end
end
if ~isempty(this.mPosData)
    if strcmp(this.plotMode, '3d')
        plot3(this.mPosData(2,:),this.mPosData(3,:),this.mPosData(4,:), 'ro','markerfacecolor',[0 0 0]);
    else
        plot(this.mPosData(1,:) * this.tStep, this.mPosData(this.plot2dDim+1,:), 'ro','markerfacecolor',[0 0 0]);
    end
end

hold off;

if(this.plotPDF && strcmp(this.plotMode, '2d'))
    subplot(2,1,2);
    
    if isa(representation, 'representations.Ensemble')
        [bw,p,x] = tools.kde( representation.ensembleM(this.plot2dDim,:));
        plot(x, p, '-');
        ylim([0 max(p)*1.1]);
        xlim([mdown mup]);
        set(gca,'YTickLabel','');
        set(gca,'YTick',[]);
    elseif isa(representation, 'representations.PCE')
        vals = pce_evaluate(representation.coefficients, representation.basis, representation.sample);
        [bw,p,x] = tools.kde(vals(this.plot2dDim,:));
        plot(x, p, '-');
        ylim([0 max(p)*1.1]);
        xlim([mdown mup]);
        set(gca,'YTickLabel','');
        set(gca,'YTick',[]);
    end
end
end