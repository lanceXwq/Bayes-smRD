function p = ucg(time_eff, range, variance)
    time_eff2 = time_eff ./ (sqrt(2 * variance));
    range_eff = range ./ (sqrt(2 * variance));
    p = 1 ./ (2 * range) .* (erf(-time_eff2 + range_eff) - erf(-time_eff2));
