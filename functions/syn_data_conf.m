function syn_data_conf(phys, traj, K, hyparam)
    N = traj.num_ptcl;
    c = zeros(K, 1, N, 'int8');
    % sample the very first conformational state
    if ~isempty(hyparam)
        c(1, 1, :) = randc(repmat(hyparam.init_conf, [1, 1, N]));
    else
        % TODO
        c(1, 1, :) = randc(repmat([0.5, 0.5], [1, 1, N]));
    end

    % sample the state at each pulse time
    p = zeros(1, phys.num_spc, N);

    for k = 2:K

        for n = 1:N
            p(1, :, n) = phys.trans_prob(c(k - 1, 1, n), :);
        end

        c(k, 1, :) = randc(p);
    end

    traj.conf = permute(c, [1, 3, 2]);
