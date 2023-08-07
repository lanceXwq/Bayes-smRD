classdef Sample

    properties
        num_roi uint32 = 1
        phys = Physics()
        traj = Trajectories()
    end

    properties (Dependent)
        log_post_full
        log_post_no_spat
        log_post_no_traj
        mem_usage
    end

    methods

        function obj = Sample(varargin)

            if nargin == 1
                obj.num_roi = varargin{1};
                obj.traj(obj.num_roi, 1) = Trajectories();
            end

        end

        function obj = calculate_log_post(obj, detn_all, param, hyparam, trans_num)
            obj.phys.calculate_log_post(hyparam);

            for idx_roi = 1:obj.num_roi
                obj.traj(idx_roi).calculate_log_post(...
                    obj.phys, detn_all(idx_roi), param, hyparam, trans_num(:, :, idx_roi));
            end

        end

        function value = get.log_post_full(obj)
            value = 0;

            for idx_roi = 1:obj.num_roi
                value = value + obj.traj(idx_roi).log_lh_obs;
                value = value + obj.traj(idx_roi).log_prior_conf;
                value = value + obj.traj(idx_roi).log_prior_spat;
                value = value + obj.traj(idx_roi).log_prior_fluor;
            end

            value = value + sum(obj.phys.log_prior_diff_coef);
            value = value + obj.phys.log_prior_exc;
            value = value + sum(obj.phys.log_prior_FRET);
        end

        function value = get.log_post_no_spat(obj)
            value = 0;

            for idx_roi = 1:obj.num_roi
                value = value + obj.traj(idx_roi).log_lh_obs;
                value = value + obj.traj(idx_roi).log_prior_conf;
                value = value + obj.traj(idx_roi).log_prior_fluor;
            end

            value = value + sum(obj.phys.log_prior_diff_coef);
            value = value + obj.phys.log_prior_exc;
            value = value + sum(obj.phys.log_prior_FRET);
        end

        function value = get.log_post_no_traj(obj)
            value = 0;

            for idx_roi = 1:obj.num_roi
                value = value + obj.traj(idx_roi).log_lh_obs;
                value = value + obj.traj(idx_roi).log_prior_fluor;
            end

            value = value + sum(obj.phys.log_prior_diff_coef);
            value = value + obj.phys.log_prior_exc;
            value = value + sum(obj.phys.log_prior_FRET);
        end

        function tot_size = get.mem_usage(obj)
            temp = obj.num_roi;
            prop_info = whos('temp');
            tot_size = prop_info.bytes;
            names = fieldnames(obj.traj);

            for idx = 1:numel(names)
                temp = [obj.traj.(names{idx})];
                prop_info = whos('temp');
                tot_size = tot_size + prop_info.bytes;
            end

            names = fieldnames(obj.phys);

            for idx = 1:numel(names)
                temp = [obj.phys.(names{idx})];
                prop_info = whos('temp');
                tot_size = tot_size + prop_info.bytes;
            end

        end

    end

end
