classdef GroundTruth < Sample

    properties
        emis Emission = Emission()
        hyparam HyperParameter = HyperParameter()
    end

    methods

        function obj = GroundTruth(num_roi)
            obj@Sample(num_roi);
            obj.emis(num_roi, 1) = Emission();
        end

        %{
        function num_emis = photonStat(obj, idx_roi)
            % Only works for single molecule.
            % All stats are done with emis, not detn.
            chan = zeros(length(obj.time_pulse{idx_roi}), 1);
            ubknd = false(length(obj.time_pulse{idx_roi}), 1);
            sig = false(length(obj.time_pulse{idx_roi}), 1);
            chan(obj.idx_emis{idx_roi}) = obj.emis(idx_roi).chan;
            sig(obj.idx_sig_pulse{idx_roi}) = true;
            ubknd(obj.idx_emis{idx_roi}(obj.idx_ubknd{idx_roi})) = true;
            num_emis = zeros(max(obj.traj.conf{idx_roi}, [], 'all'), 4);

            for idx_conf = 1:size(num_emis, 1)
                num_emis(idx_conf, 1) = nnz((chan == 1) ...
                    & (obj.traj.conf{idx_roi} == idx_conf) & sig);
                num_emis(idx_conf, 2) = nnz((chan == 2) ...
                    & (obj.traj.conf{idx_roi} == idx_conf) & sig);
                num_emis(idx_conf, 3) = nnz((chan == 1) ...
                    & (obj.traj.conf{idx_roi} == idx_conf) & ubknd);
                num_emis(idx_conf, 4) = nnz((chan == 2) ...
                    & (obj.traj.conf{idx_roi} == idx_conf) & ubknd);
            end

        end

        %}

        function value = idx_sig_donor(obj, idx_roi)
            value = obj.emis(idx_roi).src > 0 ...
                & obj.emis(idx_roi).chan == 1;
        end

        function value = idx_sig_accep(obj, idx_roi)
            value = obj.emis(idx_roi).src > 0 ...
                & obj.emis(idx_roi).chan == 2;
        end

        %{
        function value = idx_emis_donor(obj, idx_roi)
            value = obj.traj(idx_roi).idx_ubknd_donor ...
                | obj.idx_sig_donor(idx_roi);
        end

        function value = idx_emis_accep(obj, idx_roi)
            value = obj.traj(idx_roi).idx_ubknd_accep ...
                | obj.idx_sig_accep(idx_roi);
        end

        function value = idx_xtalk(obj, idx_roi)
            value = obj.emis(idx_roi).chan ~= obj.chan_detn{idx_roi};
        end

        %}
    end

end
