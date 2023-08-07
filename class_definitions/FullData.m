classdef FullData < handle
    % this is a class that holds all the informaton of the data from Ben Schuler's lab

    properties
        unit_macro = 's'
        unit_micro = 'ns'
        start double
        pulse_period double
        micro_resolution double

        time_M double
        time_m_detected double
        time_m_real double

        detector uint8
        idx_duplicate logical

        alpha = 0.042
        RCM = [1 -0.068 0 0; -0.0021 1.15365 0 0; 0 0 1.26655 -0.0679; 0 0 -0.0033 1.17389]
        gamma = 1

        rBGD double
        rBGA double
    end

    properties (Dependent)
        idx_blue
        idx_orng
        idx_donor
        idx_accep
        idx_donor_blue
        idx_accep_blue
        idx_donor_orng
        idx_accep_orng
        chan
    end

    methods

        function obj = FullData(filename)
            obj.start = 1e9 * h5readatt(...
                filename, ...
                '/InterPhotonTimes', ...
                'StartTimeInSeconds');
            obj.pulse_period = 1e-3 * double(h5readatt(...
                filename, ...
                '/InterPhotonTimes', ...
                'MacroTimeUnitInPicoSeconds'));
            obj.micro_resolution = 1e-3 * double(h5readatt(...
                filename, ...
                '/MicroTimes', ...
                'MicroTimeTimeUnitInPicoSeconds'));

            % overall detection
            pulse_interphoton = h5read(...
                filename, ...
                '/InterPhotonTimes');
            obj.time_M = 1e-9 * obj.pulse_period * cumsum(double(pulse_interphoton));
            obj.time_m_detected = obj.micro_resolution * double(h5read(...
                filename, ...
                '/MicroTimes'));
            obj.time_m_real = obj.time_m_detected ...
                - (obj.time_m_detected > obj.pulse_period / 2) .* obj.pulse_period / 2;
            obj.detector = h5read(...
                filename, ...
                '/Detectors');

            % locate all the cases where multiple photons arrive at detectors within one pulse period
            obj.idx_duplicate = pulse_interphoton(2:end) == 0;
            obj.idx_duplicate = [0; obj.idx_duplicate] | [obj.idx_duplicate; 0];
        end

        function bursts = hist_blue(obj, binWidth, varargin)

            if nargin > 2
                threshold = varargin{1};

                if nargin == 4
                    time = obj.time_M(obj.idx_blue ...
                        & obj.time_M >= varargin{2}(1) ...
                        & obj.time_M <= varargin{2}(2));

                else
                    time = obj.time_M(obj.idx_blue);
                end

            end

            figure
            h = histogram(time, 'BinWidth', binWidth);
            idx = find(h.Values >= threshold);
            bursts(:, 1) = h.BinEdges(idx);
            bursts(:, 2) = h.BinEdges(idx + 1);
            bursts(:, 3) = h.Values(idx);
            [~, idx] = sort(bursts(:, 3), 'descend');
            bursts = bursts(idx, :);
        end

        function bursts = hist_orng(obj, binWidth, varargin)

            if nargin > 2
                threshold = varargin{1};

                if nargin == 4
                    time = obj.time_M(obj.idx_orng ...
                        & obj.time_M >= varargin{2}(1) ...
                        & obj.time_M <= varargin{2}(2));

                else
                    time = obj.time_M(obj.idx_orng);
                end

            end

            figure
            h = histogram(time, 'BinWidth', binWidth);
            idx = find(h.Values >= threshold);
            bursts(:, 1) = h.BinEdges(idx);
            bursts(:, 2) = h.BinEdges(idx + 1);
            bursts(:, 3) = h.Values(idx);
            [~, idx] = sort(bursts(:, 3), 'descend');
            bursts = bursts(idx, :);
        end

        function check(obj, time_interval)
            color_d = [0, 204/255, 51/255];
            color_a = [1, 51/255, 0];

            num_fig = length(findobj('type', 'figure'));
            idx_full = [...
                        obj.idx_donor_blue, ...
                        obj.idx_accep_blue, ...
                        obj.idx_donor_orng, ...
                        obj.idx_accep_orng];

            for idx_int = size(time_interval, 1):-1:1

                schuler_correction(obj, time_interval(idx_int, :))

                idx = idx_full ...
                    & obj.time_M >= time_interval(idx_int, 1) ...
                    & obj.time_M <= time_interval(idx_int, 2);
                figure(num_fig + idx_int)
                tiledlayout(2, 1)
                ax1(idx_int) = nexttile;
                stem(obj.time_M(idx(:, 1)), obj.time_m_real(idx(:, 1)), ...
                    'Color', color_d, ...
                    'Marker', 'none'); hold on
                text(1, 0.95, num2str(nnz(idx(:, 1))), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', ...
                    'Color', color_d, ...
                    'FontSize', 14)
                text(0, 0.95, [num2str(mean(obj.time_m_real(idx(:, 1)))), 'ns'], ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'left', ...
                    'Color', color_d, ...
                    'FontSize', 14)
                stem(obj.time_M(idx(:, 2)), -obj.time_m_real(idx(:, 2)), ...
                    'Color', color_a, ...
                    'Marker', 'none'); hold off
                text(1, 0.05, num2str(nnz(idx(:, 2))), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', ...
                    'Color', color_a, ...
                    'FontSize', 14)
                text(0, 0.05, [num2str(mean(obj.time_m_real(idx(:, 2)))), 'ns'], ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'left', ...
                    'Color', color_a, ...
                    'FontSize', 14)
                title('blue laser'); ylabel('(ns)')
                ax2(idx_int) = nexttile;
                stem(obj.time_M(idx(:, 3)), obj.time_m_real(idx(:, 3)), ...
                    'Color', color_d, ...
                    'Marker', 'none'); hold on
                text(1, 0.95, num2str(nnz(idx(:, 3))), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', ...
                    'Color', color_d, ...
                    'FontSize', 14)
                text(0, 0.95, [num2str(mean(obj.time_m_real(idx(:, 3)))), 'ns'], ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'left', ...
                    'Color', color_d, ...
                    'FontSize', 14)
                stem(obj.time_M(idx(:, 4)), -obj.time_m_real(idx(:, 4)), ...
                    'Color', color_a, ...
                    'Marker', 'none'); hold off
                text(1, 0.05, num2str(nnz(idx(:, 4))), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', ...
                    'Color', color_a, ...
                    'FontSize', 14)
                text(0, 0.05, [num2str(mean(obj.time_m_real(idx(:, 4)))), 'ns'], ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'left', ...
                    'Color', color_a, ...
                    'FontSize', 14)
                xlabel(['macro time (', obj.unit_macro, ')']);
                title('orange laser'); ylabel('(ns)')
                linkaxes([ax1(idx_int), ax2(idx_int)]); xlim(time_interval(idx_int, :))
            end

        end

        function add_detection(obj, exp_data, rad_A, time_interval)

            if isempty(exp_data.param.pulse_period)
                exp_data.param.pulse_period = obj.pulse_period;
            elseif exp_data.param.pulse_period ~= obj.pulse_period
                disp('inconsistent pulse period')
                return
            end

            chan_full = obj.chan;
            exp_data.detn(exp_data.num_roi + size(time_interval, 1)) = ExperimentalDetn;

            for idx_int = 1:size(time_interval, 1)
                %! this is where I dicide how to treat the data
                % remove the detection when multiple photons hit detectors with in one pulse period
                idx_photon = find(...
                    obj.time_M >= time_interval(idx_int, 1) ...
                    & obj.time_M <= time_interval(idx_int, 2) ...
                    & obj.idx_blue ...
                    & ~obj.idx_duplicate);
                idx_pulse = uint64(obj.time_M(idx_photon) / (obj.pulse_period / 1e9));
                idx_pulse = idx_pulse - idx_pulse(1) + 1;

                n = exp_data.num_roi + idx_int;
                exp_data.detn(n).idx = idx_pulse - idx_pulse(1) +1;
                exp_data.detn(n).start_time = obj.time_M(idx_photon(1)) * 1e9;
                exp_data.detn(n).duration = obj.time_M(idx_photon(end)) * 1e9 - exp_data.detn(n).start_time;
                exp_data.detn(n).chan = chan_full(idx_photon);
                exp_data.detn(n).time_arriv = obj.time_m_real(idx_photon);

                exp_data.detn(n).prep_pulse(exp_data.param.pulse_period);
                exp_data.detn(n).prep_detn(rad_A, exp_data.param);
            end
            exp_data.num_roi = exp_data.num_roi + size(time_interval, 1);
        end

        function schuler_correction(obj, time_interval)
            nBGD = obj.rBGD' * diff(time_interval) * 1000;
            nBGA = obj.rBGA' * diff(time_interval) * 1000;

            idx = obj.time_M >= time_interval(1) ...
                & obj.time_M <= time_interval(2);
            
            nDraw = zeros(4, 1);
            for d = 1:4
                nDraw(d) = nnz((obj.detector == d) & idx & obj.idx_blue);
            end

            nAraw = zeros(4, 1);
            for d = 1:2:4
                nAraw(d) = nnz((obj.detector == d) & idx & obj.idx_orng);
            end
            
            nD = obj.RCM * (nDraw - nBGD);
            nA = nAraw - nBGA;

            nDA = nD(1) + nD(3) - obj.alpha * sum(nD);
            nDD = nD(2) + nD(4);
            nAA = nA(1) + nA(3);

            T = nDA / (nDA + nDD)
            S = (nDA + nDD)/ (nDA + nDD + obj.gamma * nAA)
        end

        function value = get.idx_blue(obj)
            value = (obj.time_m_detected / obj.pulse_period) <= 0.5;
        end

        function value = get.idx_orng(obj)
            value = (obj.time_m_detected / obj.pulse_period) > 0.5;
        end

        function value = get.idx_donor(obj)
            value = obj.detector == 2 | obj.detector == 4;
        end

        function value = get.idx_accep(obj)
            value = obj.detector == 1 | obj.detector == 3;
        end

        function value = get.idx_donor_blue(obj)
            value = obj.idx_donor & obj.idx_blue;
        end

        function value = get.idx_accep_blue(obj)
            value = obj.idx_accep & obj.idx_blue;
        end

        function value = get.idx_donor_orng(obj)
            value = obj.idx_donor & obj.idx_orng;
        end

        function value = get.idx_accep_orng(obj)
            value = obj.idx_accep & obj.idx_orng;
        end

        function value = get.chan(obj)
            value = uint8(mod(obj.detector, 2)) + 1;
        end

    end

end
