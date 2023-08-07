function P = rate2prob(Q)

    if size(Q, 1) ~= 2
        P = expm(Q);
    else
        P = zeros(2);
        t = exp(Q(1, 1) + Q(2, 2));
        P(1, 1) = (Q(1, 2) * t + Q(2, 1)) / (Q(1, 2) + Q(2, 1));
        P(1, 2) = (Q(1, 2) - Q(1, 2) * t) / (Q(1, 2) + Q(2, 1));
        P(2, 2) = 1 + t - P(1, 1);
        P(2, 1) = 1 - t - P(1, 2);
    end
