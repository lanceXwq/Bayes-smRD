function sampler_fluor(phys, traj, detn, param, hyparam)
    fluor = [true, true; true, false; false, true; false, false];
    log_post = zeros(1, 4, 'double');
    exc_prob = zeros([size(traj.exc_prob), 4], 'double');

    for idx = 1:4

        if any(fluor(idx, :) ~= traj.fluor)
            exc_prob(:, :, idx) = prob_exc(phys.exc, fluor(idx, :), traj.PSF, param.pulse_width);
        else
            exc_prob(:, :, idx) = traj.exc_prob;
        end

        log_post(idx) = sum(log_prob_obs(...
            phys.FRET, phys.rad, traj.conf, exc_prob(:, :, idx), fluor(idx, :), ...
            detn, param));
    end

    log_post = log_post + log(hyparam.fluor);
    idx = randc(exp(log_post - max(log_post)));
    traj.fluor = fluor(idx, :);
    traj.exc_prob = exc_prob(:, :, idx);
