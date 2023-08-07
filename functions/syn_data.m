function [gnd_truth, exp_data] = syn_data(gnd_truth, exp_data)

    for idx_roi = 1:gnd_truth.num_roi
        traj = gnd_truth.traj(idx_roi);
        emis = gnd_truth.emis(idx_roi);
        detn = exp_data.detn(idx_roi);
        detn.prep_pulse(exp_data.param.pulse_period);

        % generate chemical trajectory
        syn_data_conf(gnd_truth.phys, traj, detn.pulse_num, gnd_truth.hyparam);
        % generate spatial trajectory
        syn_data_spat(gnd_truth.phys, traj, exp_data.param, gnd_truth.hyparam);
        % generate emission channel
        syn_data_emis(gnd_truth.phys, traj, exp_data.param, emis);
        % generate photon arrival data
        syn_data_detn(traj, exp_data.param, emis, detn);
        
        detn.prep_detn(gnd_truth.phys.rad(2), exp_data.param);
    end

    % Data visualization
    syn_data_vis(gnd_truth, exp_data);
