function accep = sampler_FRET(phys, traj_all, exp_data, hyparam, accep)

    M = phys.num_spc;
    FRET = phys.FRET ./ hyparam.FRET_prop.' .* randg(hyparam.FRET_prop).';

    log_lh_prop = zeros(1, M, 'double');
    log_lh_old = zeros(1, M, 'double');

    for idx_roi = 1:exp_data.num_roi
        traj = traj_all(idx_roi);
        detn = exp_data.detn(idx_roi);
        log_lh_prop_temp = log_prob_obs_FRET(...
            FRET, phys.rad, traj.conf(detn.idx), traj.exc_prob, traj.fluor, ...
            detn, exp_data.param);
        log_lh_old_temp = log_prob_obs_FRET(...
            phys.FRET, phys.rad, traj.conf(detn.idx), traj.exc_prob, traj.fluor, ...
            detn, exp_data.param);

        for m = 1:M
            idx = traj.conf(detn.idx) == m;
            log_lh_prop(m) = log_lh_prop(m) + sum(log_lh_prop_temp(idx));
            log_lh_old(m) = log_lh_old(m) + sum(log_lh_old_temp(idx));
        end
    end
    
    for m = M:-1:1
    	log_r(m) = log_lh_prop(m) - log_lh_old(m) ...
            + (2 * hyparam.FRET_prop(m) - hyparam.FRET(m, 1)) * log(phys.FRET(m) / FRET(m)) ...
            + hyparam.FRET(m, 1) / hyparam.FRET(m, 2) * (phys.FRET(m) - FRET(m)) ...
            + hyparam.FRET_prop(m) * (FRET(m) / phys.FRET(m) - phys.FRET(m) / FRET(m));
        if log(rand) < log_r(m)
            accep(m, 1) = accep(m, 1) + 1;
            phys.FRET(m) = FRET(m);
        end
    end
    
    accep(:, 2) = accep(:, 2) + 1;
