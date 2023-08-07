classdef ExperimentalParam < handle

    properties
        pulse_period double
        pulse_width double
        %confocal_dim double {mustBeNonnegative} = [0.3, 0.3, 1.5]
        confocal_dim double {mustBeNonnegative} = [0.841, 0.841, 1.8]
        container_dim double {mustBeNonnegative} = [9, 9, 9]

        ubknd_rate double {mustBeNonnegative} = [0, 0]
        ubknd_prob double = [0, 0]

        exc_ratio double {mustBeNonnegative} = [1, 0]

        xtalk_prob double {mustBeNonnegative, mustBeNormalized} ...
            = [0.94, 0.06; 0.01, 0.99]

        irf_offset double {mustBeNonnegative} = 0.924
        irf_stddev double {mustBeNonnegative} = 0.224
        irf_var double {mustBeNonnegative} = 0.924^2
    end

    methods

        function prep(obj)
            obj.ubknd_prob = expm1(obj.ubknd_rate * obj.pulse_period);
            obj.irf_var = obj.irf_stddev^2;
        end

    end

end
