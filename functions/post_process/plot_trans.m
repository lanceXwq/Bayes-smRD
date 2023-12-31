function plot_trans(prob, bounds, MAP, prior_param, varargin)

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
    tiledlayout(1, 2);
    nexttile;
    h1 = histogram(prob(1, 2:2:end), ...
        'Normalization', 'pdf', ...
        'FaceColor', colors(1, :), ...
        'FaceAlpha', 1);
    hold on
    area([bounds(1, 1), bounds(1, 3)], [max(h1.Values), max(h1.Values)], ...
        'FaceColor', [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    uistack(h1, 'top')
    x = linspace(0, max(h1.BinEdges), 3001);
    plot(x, betapdf(x, prior_param(1, 2), prior_param(1, 1)), ...
        'Color', colors(2, :), ...
        'LineWidth', line_width)
    plot([MAP(1, 2), MAP(1, 2)], [0, max(h1.Values)], ...
        'Color', colors(4, :), ...
        'LineWidth', line_width)

    if nargin == 5
        plot([varargin{1}(1, 2), varargin{1}(1, 2)], [0, max(h1.Values)], ...
            'Color', colors(3, :), ...
            'LineWidth', line_width)
    end

    hold off
    ylabel('PDF')
    xlabel('Trans. prob.')
    ylim([0, max(h1.Values)])
    xlim([0, max(h1.BinEdges)])
    box off

    nexttile;
    h2 = histogram(prob(2, 1:2:end), ...
        'Normalization', 'pdf', ...
        'FaceColor', colors(1, :), ...
        'FaceAlpha', 1);
    hold on
    area([bounds(2, 1), bounds(2, 3)], [max(h2.Values), max(h2.Values)], ...
        'FaceColor', [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    uistack(h2, 'top')
    x = linspace(0, max(h2.BinEdges), 3001);
    plot(x, betapdf(x, prior_param(2, 1), prior_param(2, 2)), ...
        'Color', colors(2, :), ...
        'LineWidth', line_width)
    plot([MAP(2, 1), MAP(2, 1)], [0, max(h2.Values)], ...
        'Color', colors(4, :), ...
        'LineWidth', line_width)

    if nargin == 5
        plot([varargin{1}(2, 1), varargin{1}(2, 1)], [0, max(h1.Values)], ...
            'Color', colors(3, :), ...
            'LineWidth', line_width)
    end

    hold off
    ylabel('PDF')
    xlabel('Trans. prob.')
    ylim([0, max(h2.Values)])
    xlim([0, max(h2.BinEdges)])
    box off

    lgd = legend(labels, ...
        'NumColumns', num_columns, ...
        'FontSize', 11);
    lgd.Box = 'off';
    lgd.Layout.Tile = 'North';
