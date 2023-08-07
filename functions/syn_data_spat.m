function syn_data_spat(phys, traj, param, hyparam)
    [K, N] = size(traj.conf);

    % establish the sequence of diffusion coefficients at each pulse
    D = phys.diff_coef(traj.conf);

    % make sure the result is a column vector
    if N == 1
        D = D.';
    end

    % sample the middle position
    if ~isempty(hyparam)
        mid_pos = hyparam.init_pos_mean + sqrt(hyparam.init_pos_var) .* randn(1, 3, N);
    else
        mid_pos = zeros(1, 3, N);
    end

    idx_mid = round(K / 2);
    % calculate all the displacements between any two adjacent pulses
    displ_prev = sqrt(2 * D(idx_mid - 1:-1:1, 1, :) * param.pulse_period) ...
        .* randn(idx_mid - 1, 3, N);
    displ_post = sqrt(2 * D(idx_mid:end - 1, 1, :) * param.pulse_period) ...
        .* randn(K - idx_mid, 3, N);
    % add displacements up to form the trajectory
    x = [flipud(mid_pos - cumsum(displ_prev)); cumsum([mid_pos; displ_post])];

    %x = ones(K, 3, 'double') / 4;

    traj.update_spat(x, param.confocal_dim, false);
    traj.exc_prob = prob_exc(phys.exc, traj.fluor, traj.PSF, param.pulse_width);

    % boundary condition
    %{
    smp.roi.spat = smp.roi.spat ...
        - fix(smp.roi.spat ./ exp_data.container_dim * 2) ...
        .* exp_data.container_dim / 2;
    %}
