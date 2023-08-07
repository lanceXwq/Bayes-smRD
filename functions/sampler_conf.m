function sampler_conf(phys, traj, detn, param, hyparam)
    K = size(traj.spat, 1);
    M = phys.num_spc;
    log_lh = zeros(K, M); % log likelihood

    for m = 1:M
        % calculate the photon arrival time and channel of detection part
        log_lh(detn.idx, m) = log_prob_obs_conf( ...
        phys.FRET, phys.rad, m, traj.exc_prob(detn.idx), traj.fluor, detn, param);
        % calculate the diffusion part
        log_lh(:, m) = log_lh(:, m) ...
        + log_prob_spat(phys.diff_coef, m, traj.spat, param);
    end

    % construct the filters and sample from them
    traj.conf = forwback_mex(exp(log_lh - max(log_lh, [], 2)), hyparam.init_conf, phys.trans_prob, rand(K, 1));
