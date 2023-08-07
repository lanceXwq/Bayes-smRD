function log_prob = log_prob_obs_FRET(FRET, rad, conf, exc_prob, fluor, detn, param)
    K = length(detn.time_eff);
    esc = esc_eff(FRET, rad(1), fluor, conf);
    emr = esc - rad(2);
    prob_donor_emis = rad(1) ./ esc;
    prob_FRET = 1 - prob_donor_emis;

    ecg_esc_prob = prob_ecg(detn.time_eff, esc, param.irf_var);

    log_prob = ones(K, 1, 'double');
    log_prob(detn.idx_donor) ...
        = sum(...
        exc_prob(detn.idx_donor, 1:2:end) .* prob_donor_emis(detn.idx_donor) .* param.xtalk_prob(1, 1) ...
        .* ecg_esc_prob(detn.idx_donor), 2) ...
        + sum(...
        exc_prob(detn.idx_donor, 1:2:end) .* prob_FRET(detn.idx_donor) .* param.xtalk_prob(2, 1) ./ emr(detn.idx_donor) ...
        .* (esc(detn.idx_donor) .* detn.ecg_A_prob(detn.idx_donor) - rad(2) .* ecg_esc_prob(detn.idx_donor)), 2) ...
        + sum(...
        exc_prob(detn.idx_donor, 2:2:end) .* param.xtalk_prob(2, 1), 2) .* detn.ecg_A_prob(detn.idx_donor) ...
        + param.ubknd_prob(1) / param.pulse_period;
    log_prob(detn.idx_accep) ...
        = sum(...
        exc_prob(detn.idx_accep, 1:2:end) .* prob_donor_emis(detn.idx_accep) .* param.xtalk_prob(1, 2) ...
        .* ecg_esc_prob(detn.idx_accep), 2) ...
        + sum(...
        exc_prob(detn.idx_accep, 1:2:end) .* prob_FRET(detn.idx_accep) .* param.xtalk_prob(2, 2) ./ emr(detn.idx_accep) ...
        .* (esc(detn.idx_accep) .* detn.ecg_A_prob(detn.idx_accep) - rad(2) .* ecg_esc_prob(detn.idx_accep)), 2) ...
        + sum(...
        exc_prob(detn.idx_accep, 2:2:end) .* param.xtalk_prob(2, 2), 2) .* detn.ecg_A_prob(detn.idx_accep) ...
        + param.ubknd_prob(2) / param.pulse_period;

    log_prob = log(log_prob);
