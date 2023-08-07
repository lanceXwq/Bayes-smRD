function plot_FRET_eff(rate, bounds, MAP, prior_param, varargin)

    colors = linspecer(4);
    labels = {...
            '95% conf.', ...
            'Post. prob. dist.', ...
            'Prior dist.', ...
            'MAP'};

    if nargin == 5
        labels{end + 1} = 'Real';
    end

    % TODO input donor emission rate
    eff = rate ./ (rate + 0.244);
    eff_bounds = bounds ./ (bounds + 0.244);
    eff_MAP = MAP ./ (MAP + 0.244);

    num_columns = numel(labels);
    line_width = 2;
    figure
    tiledlayout(1, size(eff, 1));

    for m = 1:size(eff, 1)
        nexttile;
        h = histogram(eff(m, :), ...
            'Normalization', 'pdf', ...
            'FaceAlpha', 1, ...
            'FaceColor', colors(1, :));
        hold on
        area([eff_bounds(m, 1), eff_bounds(m, 3)], [max(h.Values), max(h.Values)], ...
            'FaceColor', [0, 0, 0.6], ...
            'FaceAlpha', 0.1, ...
            'EdgeColor', 'none')
        uistack(h, 'top')

        x = linspace(min(h.BinEdges), max(h.BinEdges), 3001);

        plot(x, eff_prior(x, prior_param(m, 1), prior_param(m, 2) / prior_param(m, 1), 0.244), 'Color', colors(2, :), ...
            'LineWidth', line_width)
        plot([eff_MAP(m), eff_MAP(m)], [0, max(h.Values)], ...
            'Color', colors(4, :), ...
            'LineWidth', line_width)

        if nargin == 5
            plot([varargin{1}(1), varargin{1}(1)], [0, max(h.Values)], ...
                'Color', colors(3, :), ...
                'LineWidth', line_width)
        end

        hold off
        xlabel('FRET efficiency')
        ylabel('PDF')
        ylim([0, max(h.Values)])
        box off
    end

    lgd = legend(labels, ...
        'NumColumns', num_columns, ...
        'FontSize', 11);
    lgd.Box = 'off';
    lgd.Layout.Tile = 'North';
end

function p = eff_prior(eff, k, theta, rD)
    %myFun - Description
    %
    % Syntax: p = eff_prior(eff, k, theta)
    %
    % Calculate priors of FRET efficiencies.
    p = -gammaln(k) - k * log(theta) + (k - 1) * log(eff) ...
        + k * log(rD) - (k + 1) * log(1 - eff) ...
        - eff ./ (1 - eff) * rD / theta;
    p = exp(p);

end
