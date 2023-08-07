classdef HyperParameter < handle

    properties
        init_conf double {mustBeNonnegative, mustBeNormalized(init_conf)} = [0.5, 0.5]
        init_pos_mean(:, 3) double {mustBeNumeric} = [0, 0, 0]
        init_pos_var(:, 3) double {mustBeNonnegative} = [0, 0, 0]
        diff_coef double {mustBeNonnegative}
        fluor double = [0.25, 0.25, 0.25, 0.25]
        trans_prob double {mustBeNonnegative, mustBeSquare(trans_prob)}
        FRET double {mustBeNonnegative} = [3, 1.25e-1; 3, 5e-1]
        FRET_prop double {mustBeNonnegative} = [3, 3]
        exc double {mustBeNonnegative} = [10, 1e1]
    end

    properties (SetAccess = private)
        X0 (1, 3) double {mustBeNonnegative}
        inv_var (1, 3) double {mustBeNonnegative}
    end

    methods

        function obj = set_HMC_param(obj, confocal_dim)
            obj.X0 = obj.init_pos_mean ./ confocal_dim;
            obj.inv_var = (confocal_dim.^2) ./ obj.init_pos_var;
        end

        function log_prob = log_prior_init_pos(obj, x)
            log_prob = -log(obj.init_pos_var) / 2 ...
                - (x - obj.init_pos_mean).^2 ./ (2 * obj.init_pos_var);
            log_prob = sum(log_prob);
        end

        function log_prob = log_prior_init_conf(obj, c)
            log_prob = log(obj.init_conf(c));
        end

        function log_prob = log_prior_fluor(obj, f)
            f = ~f;
            log_prob = log(obj.fluor(f(1) * 2 + f(2) + 1));
        end

        function log_prob = log_prior_FRET(obj, FRET)
            log_prob(2) = log(...
                gampdf(FRET(2), obj.FRET(2, 1), obj.FRET(2, 2) / obj.FRET(2, 1)));
            log_prob(1) = log(...
                gampdf(FRET(1), obj.FRET(1, 1), obj.FRET(1, 2) / obj.FRET(1, 1)));
        end

        function log_prob = log_prior_exc(obj, exc)
            log_prob = log(gampdf(exc, obj.exc(1), obj.exc(2) / obj.exc(1)));
        end

        function log_prob = log_prior_diff_coef(obj, diff_coef)
            log_prob(2) = log(...
                invgampdf(diff_coef(2), obj.diff_coef(2, 1), obj.diff_coef(2, 2)));
            log_prob(1) = log(...
                invgampdf(diff_coef(1), obj.diff_coef(1, 1), obj.diff_coef(1, 2)));
        end

    end

end
