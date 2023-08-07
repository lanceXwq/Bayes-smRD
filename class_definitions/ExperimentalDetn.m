classdef ExperimentalDetn < handle

    properties
        start_time double = 0
        duration double
        pulse_time double
        pulse_num uint32
        chan uint8
        time_arriv double
        time_eff double
        idx uint32
        idx_donor
        idx_accep

        ecg_A_prob double
    end

    methods

        function obj = prep_pulse(obj, pulse_period)
            obj.pulse_time = obj.start_time ...
                + pulse_period * (0:round(obj.duration / pulse_period))';
            obj.pulse_num = length(obj.pulse_time);
        end

        function prep_detn(obj, rad_A, param)
            obj.idx_donor = obj.chan == 1;
            obj.idx_accep = obj.chan == 2;
            obj.time_eff = obj.time_arriv - param.irf_offset;
            obj.ecg_A_prob = prob_ecg(obj.time_eff, rad_A, param.irf_var);
        end

    end

end
