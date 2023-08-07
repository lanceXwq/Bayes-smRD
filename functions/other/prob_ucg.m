 function p = prob_ucg(time_eff, pulse_period, irf_stddev)
	time_eff2 = time_eff ./ (sqrt(2) * irf_stddev);
	range_eff = pulse_period ./ (sqrt(2) * irf_stddev);
	p = 1 ./ (2 * pulse_period) .* (erf(-time_eff2 + range_eff) - erf(-time_eff2));