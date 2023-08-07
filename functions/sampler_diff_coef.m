function sampler_diff_coef(phys, traj, param, hyparam, trans_num)
    % pre-allocation
    K = sum(trans_num, 2).';
    M = phys.num_spc;
    dist_sqrd_sum = zeros(1, M, 'double');

    for idx_roi = 1:length(traj)
        dist_sqrd = sum(diff(traj(idx_roi).spat).^2, 2);
        % sum the displacement squared for each chemical state
        for m = 1:M
            idx = traj(idx_roi).conf == m;
            dist_sqrd_sum(m) = dist_sqrd_sum(m) ...
                + sum(dist_sqrd(idx(1:end - 1)));
        end

    end

    % sample from an inverse-gamma distribution
    alpha = 1.5 * K + hyparam.diff_coef(:, 1).';
    beta = 0.25 * dist_sqrd_sum / param.pulse_period + hyparam.diff_coef(:, 2).';
    phys.diff_coef = invgamrnd(alpha, beta);
