format compact
%clear

addpath(genpath('./functions'))
addpath(genpath('./class_definitions'))

%%`
units.time = 'ns';
units.space = '\mum';

exp_data_ready = false;
exp_data_file = './exp_data.mat';
gnd_truth_file = './gnd_truth.mat';
num_roi = 6;

% set random number seed
rng(3);

%% Synthetic data generation.
if ~exist(exp_data_file, 'file') || (~exp_data_ready)
    %% Set experiment paramenters
    exp_data = ExperimentalData(num_roi);
    gnd_truth = GroundTruth(num_roi);

    exp_data.units = units;
    exp_data.param.confocal_dim = [0.3, 0.3, 1.5];
    exp_data.param.container_dim = [1e3, 1e3, 1e3];
    exp_data.param.pulse_period = 51.37;
    
    exp_data.param.pulse_width = 1e-2;
    exp_data.param.ubknd_rate = [2e-6, 2e-6];
    exp_data.param.exc_ratio = [1, 0.05];
    exp_data.param.xtalk_prob = [.94, .06; .01, .99];
    exp_data.param.irf_offset = 1.26;
    exp_data.param.irf_stddev = 0.388;
    exp_data.param.prep();

    for i = 1:num_roi
        exp_data.detn(i).start_time = 0;
        exp_data.detn(i).duration = exp_data.param.pulse_period * 6e4;
        gnd_truth.traj(i).num_ptcl = 1;
        gnd_truth.traj(i).fluor = [true, true];
    end

    gnd_truth.phys.trans_prob = [...
        0.9999, 0.0001; ...
        0.00015, 0.99985];
    gnd_truth.phys.diff_coef = [60e-9, 30e-9];
    gnd_truth.phys.FRET = [1.2e-1, 5e-1];

    gnd_truth.phys.exc = 2 * exp_data.param.exc_ratio;
    gnd_truth.phys.rad = [2.44e-1, 2.56e-1];

    syn_data(gnd_truth, exp_data);

    save(exp_data_file, 'exp_data', '-v7.3')
    save(gnd_truth_file, 'gnd_truth', '-v7.3')
else
    load(exp_data_file)
    exp_data.param.confocal_dim = [0.3, 0.3, 1.5];
end

%%
init_smp_ready = false;
init_smp_file = './init_smp.mat';
smp_hyparam_file = './smp_hyparam.mat';

%% Initial sample generation.
if ~exist(init_smp_file, 'file') || (~init_smp_ready)
    smp_hyparam = HyperParameter();
    smp_hyparam.init_pos_var = exp_data.param.confocal_dim.^2/4;
    smp_hyparam.set_HMC_param(exp_data.param.confocal_dim);

    smp_hyparam.fluor = [0.25, 0.25, 0.25, 0.25];

    smp_hyparam.init_conf = [0.5, 0.5];
    smp_hyparam.trans_prob = [...
        3+4e3, 1; ...
        1, 3+4e3];
    smp_hyparam.diff_coef = [...
        1, 1e-9; ...
        1, 1e-9];
    smp_hyparam.FRET = [...
        1, 2.44e-1; ...
        1, 2.44e-1];
    smp_hyparam.FRET_prop = [50; 50];
    smp_hyparam.exc = [1, 15];

    init_smp = Sample(exp_data.num_roi);
    init_smp.phys.trans_prob = [...
        0.999, 0.001; ...
        0.001, 0.999];
    init_smp.phys.diff_coef = [55e-9, 47e-9];
    init_smp.phys.exc = 1 * exp_data.param.exc_ratio; % excitation rate sampler cannot be tested individually directly
    init_smp.phys.FRET = [2e-1, 2.5e-1];
    
    for idx_roi = 1:exp_data.num_roi
        init_smp.traj(idx_roi).fluor = [true, true];
        syn_data_conf(...
            init_smp.phys, init_smp.traj(idx_roi), exp_data.detn(idx_roi).pulse_num, smp_hyparam);
        syn_data_spat(...
            init_smp.phys, init_smp.traj(idx_roi), exp_data.param, smp_hyparam);
    end

    save(init_smp_file, 'init_smp', '-v7.3')
    save(smp_hyparam_file, 'smp_hyparam', '-v7.3')
else
    load(init_smp_file)
    load(smp_hyparam_file)
end

%% HMC parameters
hmc_param.spat.step_base = 0.1;
hmc_param.spat.step_num = 5;
hmc_param.exc.step_base = 0.1;
hmc_param.exc.step_num = 1;

rng('shuffle');
rnd_seed = randi(1e6) - 1;
rng(rnd_seed);

%%
num_smp = 1e2 + 1;
datetime
tic
[sample, accep_count] ...
    = sampler(exp_data, init_smp, smp_hyparam, num_smp, hmc_param, 100, 100);
duration(0, 0, toc, 'Format', 'dd:hh:mm:ss')

save('./samples.mat', 'sample', '-v7.3')

for idx_roi = 1:sample(1).num_roi
    disp(['The acceptance ratio of spatial trajectory ', num2str(idx_roi), ' is ', ...
            num2str(accep_count.spat(idx_roi, 1) / accep_count.spat(idx_roi, 2)), '.']);
end
disp(['The acceptance ratio of the donor excitation rate is ', ...
        num2str(accep_count.exc(1) / accep_count.exc(2)), '.']);
for m = 1:sample(1).phys.num_spc
    disp(['The acceptance ratio of FRET rate ', num2str(m), ' is ' ...
        num2str(accep_count.FRET(m, 1) / accep_count.FRET(m, 2)), '.']);
end
disp(['The random number seed used in this run is ', num2str(rnd_seed), '.']);
