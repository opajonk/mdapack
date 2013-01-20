function [this] = plot(this, representation, varargin)


plot(this.truth, varargin{:});
mup = 20;
mdown = -10;
axis([1, this.deterministicDimension, mdown, mup]);
axis manual;
hold on;
if(this.hasMeasurement())
    f = this.measureOp();
    pos =  tools.computeMPos(this.measurementDimension, this.deterministicDimension, this.h);
    vals = f(this.truth);
    plot(pos, vals, 'ro','MarkerFaceColor',[0 0 0]);
end

persistent initialStd;
if (isempty(initialStd))
    initialStd = representation.std();
end

xmean = representation.mean();
xstd = representation.std();
plot(xmean, 'r-', 'LineWidth', 2);
plot(abs(this.truth - xmean), 'm-');
plot(xstd, 'g-');

% plot initial std and zero line for comparison
% plot(initialStd, 'k--');
plot(zeros(1,this.deterministicDimension), 'k--');

% plot +- 1 stddev
plot(xmean + xstd, 'r--');
plot(xmean - xstd, 'r--');

% plot max/min -- shows if the truth is "covered" by the ensemble
% at all
% plot(max(x,[],2), 'r--');
% plot(min(x,[],2), 'r--');

xlabel('sectors of a latitude circle');
ylabel('some atmospheric quantity, e.g. temperature')
title(sprintf('time = %6.2f', this.time));
hold off;
end