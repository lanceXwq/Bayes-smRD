function verifyTimeDistribution(data)
    
    time_U_D = data.time_emis((data.chan_act == 1) & (data.traj_conf == 1));
    time_U_A = data.time_emis((data.chan_act == 2) & (data.traj_conf == 1));
    time_B_D = data.time_emis((data.chan_act == 1) & (data.traj_conf == 2));
    time_B_A = data.time_emis((data.chan_act == 2) & (data.traj_conf == 2));
    
    rate_esc_U = data.rate_rad_donor + data.rate_FRET_U;
    rate_esc_B = data.rate_rad_donor + data.rate_FRET_B;
    
    t_max = max(data.time_emis);
    
    figure
    tiledlayout(2, 2)
    
    nexttile;
    if ~isempty(time_U_D)
        t1 = linspace(0, max(time_U_D), 1001);
        p1 = rate_esc_U * exp(-rate_esc_U * t1);
        histogram(time_U_D, 'Normalization', 'pdf')
        hold on
        plot(t1, p1, 'LineWidth', 2)
        hold off
        title('U+D')
        xlim([0, t_max])
    end
    
    
    nexttile;
    if ~isempty(time_U_A)
        t2 = linspace(0, max(time_U_A), 1001);
        p2 = (rate_esc_U * data.rate_rad_accep) / (rate_esc_U - data.rate_rad_accep)...
            * (exp(-data.rate_rad_accep * t2) - exp(-rate_esc_U * t2));
        histogram(time_U_A, 'Normalization', 'pdf')
        hold on
        plot(t2, p2, 'LineWidth', 2)
        hold off
        title('U+A')
        xlim([0, t_max])
    end
    
    nexttile;
    if ~isempty(time_B_D)
        t3 = linspace(0, max(time_B_D), 1001);
        p3 = rate_esc_B * exp(-rate_esc_B * t3);
        histogram(time_B_D, 'Normalization', 'pdf')
        hold on
        plot(t3, p3, 'LineWidth', 2)
        hold off
        title('B+D')
        xlim([0, t_max])
    end
    
    nexttile;
    if ~isempty(time_B_A)
        t4 = linspace(0, max(time_B_A), 1001);
        p4 = (rate_esc_B * data.rate_rad_accep) / (rate_esc_B - data.rate_rad_accep)...
            * (exp(-data.rate_rad_accep * t4) - exp(-rate_esc_B * t4));
        histogram(time_B_A, 'Normalization', 'pdf')
        hold on
        plot(t4, p4, 'LineWidth', 2)
        hold off
        title('B+A')
        xlim([0, t_max])
    end