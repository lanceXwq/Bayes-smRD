function r = randc(p)
    s = size(p);
    s(2) = 1;
    rand_num = rand(s);

    if size(p, 2) == 1
        r = rand_num < p;
    else
        S = cumsum(p, 2);
        [~, r] = max(rand_num .* S(:, end) < S, [], 2);
    end
