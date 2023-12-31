function plot_diff(coef, bounds, MAP, prior_param, varargin)

    colors = linspecer(4);
    labels = {...
            '95% conf.', ...
            'Post. prob. dist.', ...
            'Prior dist.', ...
            'MAP'};

    if nargin == 5
        labels{end+1} = 'Real';
    end

    num_columns = numel(labels);
    line_width = 2;
    
    figure
    tiledlayout(1, size(coef, 2))
    for m = 1:size(coef, 2)
    nexttile;
    h = histogram(coef(:, m), ...
        'Normalization', 'pdf', ...
        'FaceAlpha', 1, ...
        'FaceColor', colors(1,:));
    hold on
    area([bounds(1, m), bounds(3, m)], [max(h.Values), max(h.Values)], ...
        'FaceColor', [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    uistack(h, 'top')
    x = linspace(0, max(h.BinEdges), 1001);
    plot(x, invgampdf(x, prior_param(m, 1), 1e9 * prior_param(m, 2)), 'Color', colors(2,:), ...
        'LineWidth', line_width)
    plot(1e9 * [MAP(m), MAP(m)], [0, max(h.Values)], ...
        'Color', colors(4,:), ...
        'LineWidth', line_width)
    if nargin == 5
        plot([varargin{1}(1), varargin{1}(1)], [0, max(h.Values)], ...
            'Color', colors(3, :), ...
            'LineWidth', line_width)
    end
    hold off
    ylim([0, max(h.Values)])
    ylabel('PDF')
    xlim([min(h.BinEdges), max(h.BinEdges)])
    xlabel('Diffusion coeff. (\mum^2/s)')
    box off
    end
    %{
    nexttile;
    h2 = histogram(coef(2:2:end), ...
        'Normalization', 'pdf', ...
        'FaceAlpha', 1, ...
        'FaceColor', colors(1,:));
    hold on
    area([bounds(2, 1), bounds(2, 3)], [max(h2.Values), max(h2.Values)], ...
        'FaceColor', [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    uistack(h2, 'top')
    x = linspace(0, max(h2.BinEdges), 1001);
    plot(x, invgampdf(x, prior_param(2, 1), 1e9 * prior_param(2, 2)), 'Color', colors(2,:), ...
        'LineWidth', line_width)
    plot(1e9 * [MAP(2), MAP(2)], [0, max(h2.Values)], ...
        'Color', colors(4,:), ...
        'LineWidth', line_width)
    if nargin == 5
        plot([varargin{1}(2), varargin{1}(2)], [0, max(h.Values)], ...
            'Color', colors(3, :), ...
            'LineWidth', line_width)
    end
    hold off
    ylim([0, max(h2.Values)])
    ylabel('PDF')
    xlim([min(h2.BinEdges), max(h2.BinEdges)])
    xlabel('Diffusion coeff. (\mum^2/s)')
    box off
    %}
    lgd = legend({...
                '95% conf.', ...
                'Post. prob. dist.', ...
                'Prior dist.', ...
                'MAP'}, ...
        'NumColumns', num_columns, ...
        'FontSize', 11);
    lgd.Box = 'off';
    lgd.Layout.Tile = 'North';
