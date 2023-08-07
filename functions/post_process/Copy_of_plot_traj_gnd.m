function plot_traj_gnd(detn, conf_bounds, exc_prob_bounds, gt, factor, varargin)
    %TODO add plot for ground truth
    pulse_time = detn.pulse_time / 1e6;
    donorColor = [0, 204/255, 51/255];
    accepColor = [1, 51/255, 0];
    colors = linspecer(4);

    figure
    tiledlayout(3, 1);
    ax1 = nexttile;
    stem(pulse_time(detn.idx(detn.chan == 1)), detn.time_arriv(detn.chan == 1), ...
        'Color', donorColor, ...
        'Marker', 'none');
    hold on
    stem(pulse_time(detn.idx(detn.chan == 2)), -detn.time_arriv(detn.chan == 2), ...
        'Color', accepColor, ...
        'Marker', 'none');
    hold off
    xticklabels([])
    ylabel({'Photon'; 'arrival'; '(ns)'})
    ylim([-30, 30])
    yticks([-30, 30])
    yticklabels({'30', '30'})
    box off
    lgd = legend(ax1, {'donor', 'acceptor'}, ...
        'NumColumns', 2, ...
        'FontSize', 11);
    lgd.Layout.Tile = 'North';
    lgd.Box = 'off';

    ax2 = nexttile;
    hold on
    tt = [pulse_time', fliplr(pulse_time')];
    CC = [double(conf_bounds(:, 1))' - 0.2, fliplr(double(conf_bounds(:, 3))' + 0.2)];
    fill(tt, CC, [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    stairs(pulse_time, gt.conf, ...
        'Color', colors(3, :), ...
        'LineWidth', 1.5)
    hold off
    xticklabels([])
    ylabel({'Conformational'; 'state'})
    ylim([0, 3])
    yticks([1, 2])
    lgd = legend(ax2, {'95% conf.', 'MAP'}, ...
        'NumColumns', 3, ...
        'FontSize', 11);
    lgd.Layout.Tile = 'North';
    lgd.Box = 'off';

    ax3 = nexttile;
    hold on
    pulse_time_comp = compress(pulse_time, factor);
    tt = [pulse_time_comp', fliplr(pulse_time_comp')];
    XX = [exc_prob_bounds(:, 1)', fliplr(exc_prob_bounds(:, 3)')];
    fill(tt, XX, [0, 0, 0.6], ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', 'none')
    plot(pulse_time, gt.exc_prob(:, 1), ...
        'Color', colors(3, :), ...
        'LineWidth', 1.5)
    hold off
    ylabel({'Prob.'; 'of exc.'})

    linkaxes([ax1 ax2 ax3], 'x')
    xlim(detn.start_time / 1e6 + [0, detn.duration] / 1e6)
    xlabel('Time (ms)')
