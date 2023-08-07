function syn_data_detn(traj, param, emis, detn)

    N = size(traj.conf, 2);
    idx = emis.src <= 2 * N;
    idx_BD = emis.src == 2 * N + 1;
    idx_BA = emis.src == 2 * N + 2;

    % introduce cross-talk
    P_xtalk = [param.xtalk_prob(1, 2); param.xtalk_prob(2, 1)];
    idx_xtalk = randc(P_xtalk(emis.chan));
    detn.chan = xor(emis.chan - 1, idx_xtalk) + 1;
    detn.chan(idx_BD) = 1;
    detn.chan(idx_BA) = 2;

    detn.idx_donor = detn.chan == 1;
    detn.idx_accep = detn.chan == 2;

    detn.time_arriv = zeros(length(emis.time), 1);
    detn.time_arriv(idx) = emis.time(idx) ...
        + param.irf_offset + randn(nnz(idx), 1) * param.irf_stddev;
    detn.time_arriv(~idx) = rand(nnz(~idx), 1) * param.pulse_period;

    detn.idx = emis.idx;

    idx_neg_arriv = detn.time_arriv < 0;

    if any(idx_neg_arriv)
        disp([num2str(nnz(idx_neg_arriv)), ' arrival time(s) are negative!'])
    end

    while any(idx_neg_arriv)
        detn.time_arriv(idx_neg_arriv) = detn.time_arriv(idx_neg_arriv) + param.pulse_period;
        detn.idx(idx_neg_arriv) = detn.idx(idx_neg_arriv) - 1;
        idx_discard = detn.idx < 1;
        detn.idx(idx_discard) = [];
        detn.time_arriv(idx_discard) = [];
        idx_neg_arriv = detn.time_arriv < 0;
    end

    idx_long_arriv = detn.time_arriv > param.pulse_period;

    if any(idx_long_arriv)
        disp([num2str(nnz(idx_long_arriv)), ' arrival time(s) are greater than the pulse period!'])
    end

    while any(idx_long_arriv)
        detn.time_arriv(idx_long_arriv) = detn.time_arriv(idx_long_arriv) - param.pulse_period;
        detn.idx(idx_long_arriv) = detn.idx(idx_long_arriv) + 1;
        idx_long_arriv = detn.time_arriv > param.pulse_period;
    end

    detn.time_eff = detn.time_arriv - param.irf_offset;
