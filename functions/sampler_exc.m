function accep = sampler_exc(phys, traj, exp_data, hyparam, hmc_param, accep)

    h = hmc_param.step_base * randg;
    m = 1;
    p = randn;
    T_i = 0.5 * p^2 / m;

    exc = phys.exc(1);
    V_roi = zeros(exp_data.num_roi, 1, 'double');
    dVdexc_roi = zeros(exp_data.num_roi, 1, 'double');
    exc_prob = cell(exp_data.num_roi, 1);

    for idx_roi = 1:exp_data.num_roi
        [V_roi(idx_roi), dVdexc_roi(idx_roi)] = prob_chan_deriv_exc(...
            phys.FRET, phys.rad, ...
            traj(idx_roi).conf, traj(idx_roi).PSF, traj(idx_roi).exc_prob, traj(idx_roi).fluor, ...
            exp_data.detn(idx_roi), exp_data.param);
    end

    V_prior = -(hyparam.exc(1) - 1) * log(phys.exc(1)) ...
        + hyparam.exc(1) / hyparam.exc(2) * phys.exc(1);
    dVdexc_prior = hyparam.exc(1) / hyparam.exc(2) - (hyparam.exc(1) - 1) / phys.exc(1);

    V_i = sum(V_roi) + V_prior;
    dVdexc = sum(dVdexc_roi) + dVdexc_prior;

    for idx_hmc = 1:hmc_param.step_num
        % half-step (apply V)
        p = p - 0.5 * h .* dVdexc;
        % whole-step (apply T)
        exc = exc + h * p / m;

        if exc < 0
            exc = -exc;
            p = -p;
        end

        for idx_roi = 1:exp_data.num_roi
            exc_prob{idx_roi} = prob_exc(...
                exc * exp_data.param.exc_ratio, traj(idx_roi).fluor, traj(idx_roi).PSF, exp_data.param.pulse_width);
            [V_roi(idx_roi), dVdexc_roi(idx_roi)] = prob_chan_deriv_exc(...
                phys.FRET, phys.rad, ...
                traj(idx_roi).conf, traj(idx_roi).PSF, exc_prob{idx_roi}, traj(idx_roi).fluor, ...
                exp_data.detn(idx_roi), exp_data.param);
        end

        dVdexc = sum(dVdexc_roi) + dVdexc_prior;
        % half-step (apply V)
        p = p - 0.5 * h .* dVdexc;
    end

    % Calculate and compare Hamiltonians.
    T_f = 0.5 * p^2 / m;
    V_f = sum(V_roi) + V_prior;

    log_r = T_i + V_i - T_f - V_f;

    if log(rand) < log_r
        accep(1) = accep(1) + 1;
        phys.exc = exc * exp_data.param.exc_ratio;

        for idx_roi = 1:exp_data.num_roi
            traj(idx_roi).exc_prob = exc_prob{idx_roi};
        end

    end

    accep(2) = accep(2) + 1;
