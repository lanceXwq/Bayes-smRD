classdef Trajectories < handle & matlab.mixin.Copyable

    properties
        fluor logical = [1, 1]
        num_ptcl uint8 = 1
        conf uint8
        spat double
        coord double
        PSF double
        exc_prob double

        conf_sparse double
        exc_prob_comp double

        % log posteriors
        log_lh_obs double = -Inf
        log_prior_conf double = -Inf
        log_prior_spat double = -Inf
        log_prior_fluor double = -Inf
    end

    methods

        function update_spat(obj, x, confocal_dim, normalized)

            if normalized
                obj.spat = x .* confocal_dim;
                r = permute(sum(x.^2, 2), [1, 3, 2]);
            else
                obj.spat = x;
                r = permute(sum((x ./ confocal_dim).^2, 2), [1, 3, 2]);
            end

            obj.PSF = psf(r);
        end

        function update_spat_sigmoid(obj, x, confocal_dim, normalized, a, b)

            if normalized
                obj.spat = x .* confocal_dim;
            else
                obj.spat = x;
            end
            
            obj.PSF = prod(...
                erf(a*b*(1 + x ./ confocal_dim)) + erf(a*b*(1 - x ./ confocal_dim)), ...
                2) / 8;
        end

        function calculate_log_post(obj, phys, detn, param, hyparam, trans_num)
            obj.log_lh_obs = sum(log_prob_obs(...
                phys.FRET, phys.rad, obj.conf, obj.exc_prob, obj.fluor, detn, param));

            obj.log_prior_spat = sum(log_prob_spat(phys.diff_coef, obj.conf, obj.spat, param)) ...
                + hyparam.log_prior_init_pos(obj.spat(1, :, :));

            obj.log_prior_conf = sum(double(trans_num) .* log(phys.trans_prob), 'all') ...
                + hyparam.log_prior_init_conf(obj.conf(1, :));

            obj.log_prior_fluor = hyparam.log_prior_fluor(obj.fluor);
        end

        function value = trans_num(obj, num_spc)
            value = zeros(num_spc, 'uint32');

            trans = obj.conf(1:end - 1) + (obj.conf(2:end) - 1) * num_spc;
            c = 1;

            % count transitions between states
            for idx_col = 1:num_spc

                for idx_row = 1:num_spc
                    value(idx_row, idx_col) = nnz(c == trans);
                    c = c + 1;
                end

            end

        end

        function sparser(obj)
            obj.conf_sparse = [double(obj.conf(1)); sparse(diff(double(obj.conf)))];
        end

        function conf_full = desparser(obj)
            conf_full = uint8(cumsum(full(obj.conf_sparse)));
        end

        %{
        function value = idx_pulse_sig(obj)
            value = obj.idx_pulse_emis(obj.idx_emis_sig);
        end

        function value = idx_emis_sig(obj)
            value = obj.src_emis > 0;
        end

        function value = idx_exc_D(obj)
            value = rem(obj.src_emis, 2) == 1;
        end

        function value = idx_exc_A(obj)
            value = ~rem(obj.src_emis, 2) & obj.src_emis > 0;
        end

        function value = idx_ubknd_D(obj)
            value = obj.src_emis == -1;
        end

        function value = idx_ubknd_A(obj)
            value = obj.src_emis == -2;
        end

        function value = idx_ubknd(obj)
            value = obj.idx_ubknd_D | obj.idx_ubknd_A;
        end

        function value = num_emis(obj)
            value = length(obj.idx_pulse_emis);
        end

        function value = num_sig(obj)
            value = nnz(obj.idx_sig_emis);
        end

        function value = num_ubknd(obj)
            value = nnz(obj.idx_ubknd);
        end

        %}
    end

end
