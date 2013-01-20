function [this] = plot(this, representation, varargin)

%% plot truth
plot(this.truth, varargin{:});
            % mup = ceil(max(this.state) * 1.4);
            % mdown = floor(min(this.state) * 1.4);
            mup = 15;
            mdown = -3;
            axis([0, this.deterministicDimension, mdown, mup]);
            axis manual;
            hold on;

            %% plot measurement pos
            if(this.hasMeasurement())
                h = this.measureOp();
                pos =  h((1:this.deterministicDimension)');
                vals = h(this.truth);
                % vals = zeros(length(pos),1);
                plot(pos, vals, 'ro','MarkerFaceColor',[0 0 0]);
            end
            


%% plot ensemble
% persistent initialStd;
% if (isempty(initialStd))
%     initialStd = representation.std();
% end

xmean = representation.mean();
xstd = representation.std();

% xsmry = representation.summary();

% xmedian = xsmry(:,3);
xstd = representation.std();
plot(xmean, 'r-', 'LineWidth', 2);
plot(abs(this.truth - xmean), 'm-');
plot(xstd, 'g-');

% plot initial std and zero line for comparison
% plot(initialStd, 'k--');
plot(zeros(1,this.deterministicDimension), 'k--');

% plot percentiles
% plot(xmean-xstd, 'r:');
plot(xmean-xstd, 'r--');
plot(xmean+xstd, 'r--');
% plot(xsmry(:,5), 'r:');

% plot max/min -- shows if the truth is "covered" by the ensemble
% at all
% plot(max(x,[],2), 'r--');
% plot(min(x,[],2), 'r--');

xlabel('spatial dimension (circular)');
            ylabel('some quantity');
            title(sprintf('time = %6.2f', this.time));
            hold off;
end