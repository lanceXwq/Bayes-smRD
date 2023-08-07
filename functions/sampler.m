function [smp, accep] = sampler(exp_data, init_smp, hyparam, num_smp, hmc_param, stride, comp)
    %% Preallocation
    smp((num_smp - 1) / stride + 1, 1) = Sample(exp_data.num_roi);
    smp(1) = init_smp;

    param = exp_data.param;
    detn = exp_data.detn;
    phys = copy(init_smp.phys);
    traj = copy(init_smp.traj);

    accep_spat = zeros(exp_data.num_roi, 2, 'double');
    accep_exc = zeros(1, 2, 'double');
    accep_FRET = zeros(phys.num_spc, 2, 'double');
    trans_num = zeros(phys.num_spc, phys.num_spc, exp_data.num_roi, 'uint32');

    for idx_smp = 2:num_smp
        %tic
        for idx_roi = 1:exp_data.num_roi
            sampler_fluor(phys, traj(idx_roi), detn(idx_roi), param, hyparam);

            accep_spat(idx_roi, :) = sampler_spat(...
                phys, traj(idx_roi), detn(idx_roi), param, ...
                hyparam, hmc_param.spat, accep_spat(idx_roi, :));

            %accep_spat(idx_roi, :) = sampler_spat_sigmoid(...
            %    phys, traj(idx_roi), detn(idx_roi), param, ...
            %    hyparam, hmc_param.spat, accep_spat(idx_roi, :));

            sampler_conf(phys, traj(idx_roi), detn(idx_roi), param, hyparam);

            trans_num(:, :, idx_roi) = traj(idx_roi).trans_num(phys.num_spc);
        end

        trans_num_tot = sum(trans_num, 3);

        sampler_trans_prob(phys, trans_num_tot, hyparam);

        accep_FRET = sampler_FRET(phys, traj, exp_data, hyparam, accep_FRET);

        accep_exc = sampler_exc(phys, traj, exp_data, hyparam, hmc_param.exc, accep_exc);

        sampler_diff_coef(phys, traj, param, hyparam, trans_num_tot);
        %toc
        % Save the sample.
        if rem(idx_smp - 1, stride) == 0
            current = (idx_smp - 1) / stride + 1;
            smp(current).num_roi = exp_data.num_roi;
            smp(current).traj = copy(traj);
            smp(current).phys = copy(phys);
            smp(current).calculate_log_post(detn, param, hyparam, trans_num);
            
            for idx_roi = 1:smp(current).num_roi
                smp(current).traj(idx_roi).sparser();
                smp(current).traj(idx_roi).conf = [];
                smp(current).traj(idx_roi).spat = [];
                smp(current).traj(idx_roi).coord = [];
                smp(current).traj(idx_roi).PSF = [];
                smp(current).traj(idx_roi).exc_prob_comp ...
                    = compress(smp(current).traj(idx_roi).exc_prob(:, 1), comp);
                smp(current).traj(idx_roi).exc_prob = [];
            end

        end

    end

    accep.spat = accep_spat;
    accep.exc = accep_exc;
    accep.FRET = accep_FRET;

    smp(end).traj = traj;
    smp(end).phys = phys;
