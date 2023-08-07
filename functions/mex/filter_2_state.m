function c = filter_2_state(ll, icp, tp)
    % ll: log likelihood
    % icp: initial conformational state probability
    % tp: transition probability
    % c: conformational state trajectory

    % Construct the forward filters.
    K = size(ll, 1);
    A = zeros(K, 2, 'double');
    A(1, :) = exp(ll(1, :) + log(icp));
    A(1, :) = A(1, :) / sum(A(1, :));
    % Construct the other forward filters.
    for k = 2:K
        A(k, :) = exp(ll(k, :) + log(A(k - 1, :) * tp));
        A(k, :) = A(k, :) / sum(A(k, :));
    end

    c = zeros(K, 1, 'int8');
    % Sample the last chemical state.
    c(K) = (rand > A(K, 1)) + 1;
    % Sample the rest of the photo-states by marching backwards.
    for k = K - 1:-1:1
        p = A(k, :) .* tp(:, c(k + 1)).';
        p = p ./ sum(p);
        c(k) = (rand > p(1)) + 1;
    end
