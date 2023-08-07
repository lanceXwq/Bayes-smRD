function [V, dVdexc] = prob_chan_deriv_exc(FRET, rad, conf, PSF, exc_prob, fluor, detn, param)
    Q = 1 ./ (sum(exc_prob, 2) + sum(param.ubknd_prob) + 1);

    esc = esc_eff(FRET, rad(1), fluor, conf(detn.idx));
    prob_donor_emis = rad(1) ./ esc;
    prob_FRET = 1 - prob_donor_emis;

    ipD = detn.idx(detn.idx_donor); % idx_pulse_donor
    ipA = detn.idx(detn.idx_accep); % idx_pulse_accep

    % single particle
    % fluor type times excitation ratio times exponential
    frte = fluor .* param.exc_ratio .* (exc_prob + 1);

    q1 = exc_prob(ipD, 1) .* prob_donor_emis(detn.idx_donor) .* param.xtalk_prob(1, 1) ...
        + exc_prob(ipD, 1) .* prob_FRET(detn.idx_donor) .* param.xtalk_prob(2, 1) ...
        + exc_prob(ipD, 2) .* param.xtalk_prob(2, 1) ...
        + param.ubknd_prob(1);
    q2 = exc_prob(ipA, 1) .* prob_donor_emis(detn.idx_accep) .* param.xtalk_prob(1, 2) ...
        + exc_prob(ipA, 1) .* prob_FRET(detn.idx_accep) .* param.xtalk_prob(2, 2) ...
        + exc_prob(ipA, 2) .* param.xtalk_prob(2, 2) ...
        + param.ubknd_prob(2);

    V = -sum(log(Q)) - sum(log(q1)) - sum(log(q2));

    dVdexc = -Q .* sum(frte, 2);
    dq1dexc = frte(detn.idx_donor, 1) .* (...
        prob_donor_emis(detn.idx_donor) .* param.xtalk_prob(1, 1) ...
        + prob_FRET(detn.idx_donor) .* param.xtalk_prob(2, 1)) ...
        + frte(detn.idx_donor, 2) .* param.xtalk_prob(2, 1);
    dq2dexc = frte(detn.idx_accep, 1) .* (...
        prob_donor_emis(detn.idx_accep) .* param.xtalk_prob(1, 2) ...
        + prob_FRET(detn.idx_accep) .* param.xtalk_prob(2, 2)) ...
        + frte(detn.idx_accep, 2) .* param.xtalk_prob(2, 2);

    dVdexc(ipD) = dVdexc(ipD) + dq1dexc ./ q1;
    dVdexc(ipA) = dVdexc(ipA) + dq2dexc ./ q2;
    dVdexc = -param.pulse_width * sum(PSF .* dVdexc);
