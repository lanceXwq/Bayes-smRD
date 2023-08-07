function syn_data_vis(gnd_truth, syn_data)
    cmap = linspecer(gnd_truth.num_roi);

    for idx_roi = 1:gnd_truth.num_roi
        traj = gnd_truth.traj(idx_roi);
        emis = gnd_truth.emis(idx_roi);
        detn = syn_data.detn(idx_roi);

        pulse_time_ms = detn.pulse_time * 1e-6;
        start_time_ms = detn.start_time * 1e-6;
        duration_ms = detn.duration * 1e-6;

        figure
        tiledlayout(3, 1)

        ax1 = nexttile;
        hold on
        N = double(traj.num_ptcl);

        if mod(N, 2) == 0
            shift = (-N + 1:2:N - 1) / 10;
        else
            shift = (-(N - 1):2:(N - 1)) / 10;
        end

        for n = 1:N
            stairs(...
                pulse_time_ms, double(traj.conf(:, n)) + shift(n), ...
                'Color', cmap(idx_roi, :) ...
                );
        end

        hold off
        set(gca, 'ytick', [1, 2])
        xlim([0, duration_ms])
        ylim([0, gnd_truth.phys.num_spc + 1])
        ylabel('Conformational state')

        ax2 = nexttile;
        [x, y] = meshgrid([start_time_ms, duration_ms], linspace(0, 4, 101));
        z = exp(-2 * y.^2);
        contourf(x, y, z, 101, 'EdgeColor', 'none')
        shading interp
        colormap(flipud(gray))
        caxis([0, 2])
        view(2)
        hold on

        for n = 1:N
            plot(pulse_time_ms, -log(traj.PSF(:, n))/2, 'Color', cmap(idx_roi, :))
        end

        hold off
        title('Normalized distance')
        xlim([0, duration_ms])
        %ylim([0, 2*max(traj.exc_prob, [], 'all')])

        ax3 = nexttile;
        color_d = [0, 204/255, 51/255];
        color_a = [1, 51/255, 0];
        idx = detn.chan == [1, 2];
        T_D = pulse_time_ms(detn.idx(idx(:, 1)));
        T_A = pulse_time_ms(detn.idx(idx(:, 2)));
        t_D = detn.time_arriv(idx(:, 1));
        t_A = -detn.time_arriv(idx(:, 2));
        hold on
        stem(T_D, t_D, ...
            'Color', color_d, ...
            'Marker', 'none');
        text(1, 0.95, num2str(nnz(idx(:, 1))), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'Color', color_d, ...
            'FontSize', 14)
        text(0, 0.95, [num2str(mean(t_D)), 'ns'], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'Color', color_d, ...
            'FontSize', 14)
        stem(T_A, t_A, ...
            'Color', color_a, ...
            'Marker', 'none');
        text(1, 0.05, num2str(nnz(idx(:, 2))), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'Color', color_a, ...
            'FontSize', 14)
        text(0, 0.05, [num2str(-mean(t_A)), 'ns'], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'Color', color_a, ...
            'FontSize', 14)
        idx_xtalk = emis.chan ~= detn.chan;
        scatter(pulse_time_ms(detn.idx(idx_xtalk)), zeros(nnz(idx_xtalk), 1), 'k')
        idx_ubknd = emis.src > 2 * N; % may produce error when there is neg arriv time
        scatter(pulse_time_ms(detn.idx(idx_ubknd)), zeros(nnz(idx_ubknd), 1), 'bs')
        hold off
        xlim([0, duration_ms])
        yticks(-50:10:50);
        yticklabels([50:-10:0, 10:10:50])
        ylabel('Arrival time (ns)')
        xlabel('Pulse time (ms)')
        lgd = legend('Donor', 'Acceptor', 'Cross-talk', 'Background', 'NumColumns', 4);
        lgd.Box = 'off';
        lgd.Layout.Tile = 'north';

        linkaxes([ax1 ax2 ax3], 'x')
    end

end
