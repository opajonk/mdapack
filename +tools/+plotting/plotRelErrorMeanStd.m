function [h] = plotRelErrorMeanStd(filename, color)
    load(filename);
    xlabel('time');
    ylabel('relative error');
    h = semilogy(results.f_relError.mean, 'LineStyle', '-', 'Color', color, 'LineWidth', 2);
    hold on;
    semilogy(results.f_relError.mean + results.f_relError.std, 'LineStyle', ':', 'Color', color);
    semilogy(results.f_relError.mean - results.f_relError.std, 'LineStyle', ':', 'Color', color);
end