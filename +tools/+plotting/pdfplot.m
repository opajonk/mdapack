function [ output_args ] = pdfplot(fileName, figCap, x, varargin )

fontSize = 14;
lineWidth = 0.5;
picdims = [500 600];
figNum = 1;

fig=figure(figNum);

ha = tools.tight_subplot(1,1,[.06 .03],[.06 .01],[.1 .06]);

plot(ha(1), x,varargin{:});


%% ================== export the figure to EPS/LaTeX using modlaprint =====
% set(gca,'fontSize',fontSize);
set(fig,'Position',[0 0, picdims]);
% pos = get(fig,'position');
% set(fig,'PaperPositionMode','Auto');

modlaprint(figNum,'tex',fileName,'texoutput','texprint', ...
    'caption', figCap, ...
    'makeps',1,...
    'width',14,...
    'keepticklabels','on');
% close(fig);
end

