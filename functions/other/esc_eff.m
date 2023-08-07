function value = esc_eff(FRET, rad_donor, fluor, conf)
    value = (rad_donor + fluor(:, 2).' .* FRET(conf)).';
end