function sampler_trans_prob(phys, trans_num, hyparam)
    param = trans_num + hyparam.trans_prob;
    phys.trans_prob = randd(param);
