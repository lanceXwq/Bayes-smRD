function plot_trans_in1_2(prob, prior_param, varargin)

    colors = linspecer(4);
    labels = {...
            %'95% conf.', ...
            'Post. prob. dist.', ...
            'Prior dist.'};

    if nargin == 5
        labels{end + 1} = 'Real';
    end

    line_width = 2;
    self_prob(:, 1) = prob(1:3:end, 1);
    self_prob(:, 2) = prob(2:3:end, 2);
    %self_prob(:, 3) = prob(3:3:end, 3);

    figure
    h1 = histogram(self_prob, ...
        'Normalization', 'pdf', ...
        'FaceColor', colors(1, :), ...
        'FaceAlpha', 1);
    hold on
    x = linspace(min(h1.BinEdges), max(h1.BinEdges), 3001);
    plot(x, betapdf(x, prior_param(1, 2), prior_param(1, 1)), ...
        'Color', colors(2, :), ...
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
    xlim([0.996, 1])
    box off
