function syn_data_emis(phys, traj, param, emis)

    syn_data_emis_src(traj, param, emis);

    [K, N] = size(traj.conf);
    num_emis = length(emis.src);
    idx_exc_D = mod(emis.src, 2) & emis.src < 2 * N;
    idx_B = emis.src > 2 * N;
    idx_BA = emis.src == 2 * N + 2;
    num_exc_D = nnz(idx_exc_D);

    idx_molecule = (emis.src(idx_exc_D) + 1) / 2;
    % escape rate at each signal emission
    idx_conf = sub2ind([K, N], emis.idx(idx_exc_D), idx_molecule);
    fA = traj.fluor(idx_molecule, 2);
    esc = zeros(num_exc_D, 1, 'double') + phys.rad(1);
    esc = esc + fA .* phys.FRET(traj.conf(idx_conf)).';

    emis.time = zeros(num_emis, 1, 'double');
    emis.time(idx_exc_D) = exprnd(1 ./ esc);

    % decide if FRET happens
    idx_FRET = ~randc(phys.rad(1) ./ esc);
    idx_emis_A = false(num_emis, 1);
    idx_emis_A(idx_FRET) = true;
    idx_emis_A = idx_emis_A | (~mod(emis.src, 2) & emis.src <= 2 * N);
    emis.chan = ones(num_emis, 1, 'uint8');
    emis.chan(idx_emis_A) = 2;
    emis.chan(idx_BA) = 2;

    emis.time(idx_emis_A) = emis.time(idx_emis_A) ...
        + exprnd(1 / phys.rad(2), nnz(idx_emis_A), 1);

    emis.time(idx_B) = param.pulse_period * rand(nnz(idx_B), 1);
end

function syn_data_emis_src(traj, param, emis)
    [K, N] = size(traj.exc_prob);
    prob_src_emis = ones(K, N + 3, 'double');
    % only donor dye excited
    prob_src_emis(:, 2:2:end - 3) = traj.exc_prob(:, 1:2:end);
    % only donor dye excited
    prob_src_emis(:, 3:2:end - 2) = traj.exc_prob(:, 2:2:end);
    % only one uniform photon in either channel
    prob_src_emis(:, end - 1) = param.ubknd_prob(1);
    prob_src_emis(:, end) = param.ubknd_prob(2);

    src_emis = randc(prob_src_emis) - 1;
    emis.idx = find(src_emis ~= 0);
    emis.src = src_emis(emis.idx);
end
