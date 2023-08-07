function prob = prob_exc(exc, fluor, PSF, pulse_width)
    [K, N] = size(PSF);
    prob = zeros(K, 2 * N, 'double');
    exc_eff = pulse_width * (fluor .* exc).';
    % only donor dye excited
    prob(:, 1:2:end) = expm1(exc_eff(1, :) .* PSF);
    % only donor dye excited
    prob(:, 2:2:end) = expm1(exc_eff(2, :) .* PSF);
