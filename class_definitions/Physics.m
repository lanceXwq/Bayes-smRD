classdef Physics < handle & matlab.mixin.Copyable

    properties
        trans_prob double {mustBeNonnegative, mustBeSquare(trans_prob)}
        diff_coef double {mustBeNonnegative} = [55e-9, 47e-9]

        exc double {mustBeNonnegative} = [5, 0]
        rad double {mustBeNonnegative} = [2.44e-1, 2.56e-1]
        FRET double {mustBeNonnegative} = [1.2e-1, 5e-1]

        log_prior_diff_coef double = -Inf
        log_prior_exc double = -Inf
        log_prior_FRET double = -Inf
    end
    
    properties (Dependent)
        num_spc
    end

    methods
        
        function value = get.num_spc(obj)
            value = numel(obj.diff_coef);
        end

        function obj = prep_exc(obj, pulse_width)
            obj.exc_avg = obj.exc * pulse_width;
        end

        function value = esc(obj, c)
            value = obj.rad(1) + obj.FRET(c);

            if size(c, 2) == 1
                value = value.';
            end

        end

        function calculate_log_post(obj, hyparam)
            obj.log_prior_diff_coef = hyparam.log_prior_diff_coef(obj.diff_coef);
            obj.log_prior_exc = hyparam.log_prior_exc(obj.exc(1));
            obj.log_prior_FRET = hyparam.log_prior_FRET(obj.FRET);
        end

    end

end
