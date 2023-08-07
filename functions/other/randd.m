function r = randd(dirichlet_param)
    r = randg(dirichlet_param);
    r = r ./ sum(r, 2);