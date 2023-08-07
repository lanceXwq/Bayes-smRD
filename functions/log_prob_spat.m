function log_prob = log_prob_spat(diff_coef, conf, spat, param)

    if length(conf) == 1
        sig_sqrd = 2 * diff_coef(conf) * param.pulse_period;
    else
        sig_sqrd = 2 * diff_coef(conf(1:end - 1)).' * param.pulse_period;
    end

    log_prob = [-sum(diff(spat).^2, 2) ./ (2 * sig_sqrd) - 1.5 * log(sig_sqrd); 0];
