function plot_FRET(rate, bounds, MAP, prior_param, varargin)

    colors = linspecer(4);
    labels = {...
            '95% conf.', ...
            'Post. prob. dist.', ...
            'Prior dist.', ...
            'MAP'};

    if nargin == 5
        labels{end + 1} = 'Real';
    end

    num_columns = numel(labels);
    line_width = 2;
    figure
    tiledlayout(1, size(rate, 2));

    for m = 1:size(rate, 2)
        nexttile;
        h = histogram(rate(:, m), ...
            'Normalization', 'pdf', ...
            'FaceAlpha', 1, ...
            'FaceColor', colors(1, :));
        hold on
        area([bounds(m, 1), bounds(m, 3)], [max(h.Values), max(h.Values)], ...
            'FaceColor', [0, 0, 0.6], ...
            'FaceAlpha', 0.1, ...
            'EdgeColor', 'none')
        uistack(h, 'top')
        x = linspace(min(h.BinEdges), max(h.BinEdges), 3001);
        plot(x, gampdf(x, prior_param(m, 1), prior_param(m, 2) / prior_param(m, 1)), 'Color', colors(2, :), ...
            'LineWidth', line_width)
        plot([MAP(m), MAP(m)], [0, max(h.Values)], ...
            'Color', colors(4, :), ...
            'LineWidth', line_width)

        if nargin == 5
            plot([varargin{1}(1), varargin{1}(1)], [0, max(h.Values)], ...
                'Color', colors(3, :), ...
                'LineWidth', line_width)
        end

        hold off
        xlabel('FRET rate (ns^{-1})')
        ylabel('PDF')
        ylim([0, max(h.Values)])
        box off
    end

    %{
    nexttile;
    h2 = histogram(rate(2, :), ...
        'Normalization', 'pdf', ...
        'FaceAlpha', 1, ...
        'FaceColor', colors(1, :));
    hold on
    area([bounds(2, 1), bounds(2, 3)], [max(h2.Values), max(h2.Values)], ...
        'FaceColor', [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    uistack(h2, 'top')
    x = linspace(min(h2.BinEdges), max(h2.BinEdges), 3001);
    plot(x, gampdf(x, prior_param(2, 1), prior_param(2, 2) / prior_param(2, 1)), 'Color', colors(2, :), ...
        'LineWidth', line_width)
    plot([MAP(2), MAP(2)], [0, max(h2.Values)], ...
        'Color', colors(4, :), ...
        'LineWidth', line_width)

    if nargin == 5
        plot([varargin{1}(2), varargin{1}(2)], [0, max(h2.Values)], ...
            'Color', colors(3, :), ...
            'LineWidth', line_width)
    end

    hold off
    xlabel('low FRET rate (ns^{-1})')
    ylabel('PDF')
    ylim([0, max(h2.Values)])
    box off
    %}
    lgd = legend(labels, ...
        'NumColumns', num_columns, ...
        'FontSize', 11);
    lgd.Box = 'off';
    lgd.Layout.Tile = 'North';
