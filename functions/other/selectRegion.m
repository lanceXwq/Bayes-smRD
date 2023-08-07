function selectRegion(old_gnd, old_exp, old_init, ti, tf, suffix)
    gnd_truth = old_gnd;
    exp_data = old_exp;
    init_smp = old_init;
    
    idx_init = find(old_exp.pulse_time >= ti, 1);
    idx_final = find(old_exp.pulse_time > tf, 1) - 1;
    exp_data.start_time = old_exp.pulse_time(idx_init);
    exp_data.duration = old_exp.pulse_time(idx_final) - exp_data.start_time;
    idx = old_exp.pulse_time(exp_data.idx_detn) >= ti & old_exp.pulse_time(exp_data.idx_detn) <= tf;
    exp_data.idx_detn = exp_data.idx_detn(idx) - idx_init + 1;
    exp_data.time_detn = old_exp.time_detn(idx);
    exp_data.chan_detn = old_exp.chan_detn(idx);
    
    gnd_truth.time_pulse = old_gnd.time_pulse(idx_init:idx_final);
    gnd_truth.traj_conf = old_gnd.traj_conf(idx_init:idx_final);
    gnd_truth.traj_spat = old_gnd.traj_spat(idx_init:idx_final, :, :);
    gnd_truth.traj_norm_dist = old_gnd.traj_norm_dist(idx_init:idx_final);
    
    gnd_truth.time_emis = old_gnd.time_emis(idx);
    gnd_truth.time_detn = old_gnd.time_detn(idx);
    gnd_truth.chan_emis = old_gnd.chan_emis(idx);
    gnd_truth.chan_detn = old_gnd.chan_detn(idx);
    gnd_truth.emis_type = old_gnd.emis_type(idx);
    gnd_truth.emis_ptcl = old_gnd.emis_ptcl(idx);
    
    gnd_truth.idx_emis = old_gnd.idx_emis(idx) - idx_init + 1;
    gnd_truth.idx_sig_donor = old_gnd.idx_sig_donor(idx);
    gnd_truth.idx_sig_accep = old_gnd.idx_sig_accep(idx);
    gnd_truth.idx_xtalk = old_gnd.idx_xtalk(idx);
    
    init_smp.idx_emis = old_init.idx_emis(idx) - idx_init + 1;
    init_smp.emis_type = old_init.emis_type(idx);
    init_smp.traj_conf = old_init.traj_conf(idx_init:idx_final);
    init_smp.traj_spat = old_init.traj_spat(idx_init:idx_final, :, :);
    init_smp.traj_norm_dist = old_init.traj_norm_dist(idx_init:idx_final);
    
    save(append('syn_data_', suffix, '.mat'), 'exp_data', '-v7.3')
    save(append('gnd_truth_', suffix, '.mat'), 'gnd_truth', '-v7.3')
    save(append('init_smp_', suffix, '.mat'), 'init_smp', '-v7.3')
    