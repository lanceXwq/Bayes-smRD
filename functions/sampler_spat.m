function accep = sampler_spat(phys, traj, detn, param, hyparam, hmc_param, accep)
    phi_p1 = ((param.confocal_dim.^2).' ./ (2 * phys.diff_coef * param.pulse_period))';
    % row for species, coln for coord
    phi_p1 = phi_p1(traj.conf, :);

    X_i = traj.spat ./ param.confocal_dim;

    K = detn.pulse_num;
    h = hmc_param.step_base .* randg;
    M = ones(K, 1);
    phi_p2 = 0.5 * h .* [...
                        phi_p1(1, :) + hyparam.inv_var; ...
                        phi_p1(1:end - 2, :) + phi_p1(2:end - 1, :); ...
                        phi_p1(end - 1, :) ...
                        ];
    Mh = 2 * M ./ h;
    phi = Mh + phi_p2;
    phi_m = Mh - phi_p2;
    xi = -0.5 * h * phi_p1;

    [V_i, dVdPSF] = prob_chan_deriv_spat(...
        phys.exc, phys.FRET, phys.rad, ...
        traj.conf, traj.exc_prob, traj.fluor, ...
        detn, param);

    P = randn(K, 3);
    T_i = 0.5 * sum(P.^2 ./ M, 'all');
    U_i = V_i + 0.5 * sum((X_i(1, :) - hyparam.X0).^2 .* hyparam.inv_var) ...
        + 0.5 * sum(diff(X_i).^2 .* phi_p1(1:end - 1, :), 'all');
    dVdX = dVdPSF .* (-4 * X_i .* traj.PSF);
    X = X_i;

    for idx_hmc = 1:hmc_param.step_num
        % half-step (apply only Vq)
        P = P - 0.5 * h .* dVdX;
        % whole-step (apply M and Lq)
        X_nxt = tri_solver_3d_mex(K, phi, xi, ...
            prep_f_mex(K, X, P, h, phi_m, xi, hyparam.X0, hyparam.inv_var));
        P = Mh .* (X_nxt - X) - P;
        X = X_nxt;
        % half-step (apply only Vq)
        PSF = psf(sum(X.^2, 2));
        exc_prob = prob_exc(phys.exc, traj.fluor, PSF, param.pulse_width);
        [V_f, dVdPSF] = prob_chan_deriv_spat(...
            phys.exc, phys.FRET, phys.rad, ...
            traj.conf, exc_prob, traj.fluor, ...
            detn, param);
        dVdX = dVdPSF .* (-4 * X .* PSF);
        P = P - 0.5 * h .* dVdX;
    end

    % Calculate and compare Hamiltonians.
    T_f = 0.5 * sum(P.^2 ./ M, 'all');
    U_f = V_f + 0.5 * sum((X(1, :) - hyparam.X0).^2 .* hyparam.inv_var) ...
        + 0.5 * sum(diff(X).^2 .* phi_p1(1:end - 1, :), 'all');
    log_r = T_i + U_i - T_f - U_f;

    if log(rand) < log_r
        % accept
        accep(1) = accep(1) + 1;
        traj.update_spat(X, param.confocal_dim, true);
        traj.exc_prob = exc_prob;
    end

    accep(2) = accep(2) + 1;
