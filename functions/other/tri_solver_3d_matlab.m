function X = tri_solver_3d_matlab(K, phi, xi, f)

    X = zeros(K, 3, 'double');
    xip = zeros(K, 3, 'double');

    w = phi(1, :);
    X(1, :) = f(1, :) ./ w;

    for k = 2:K
        xip(k - 1, :) = xi(k - 1, :) ./ w;
        w = phi(k, :) - xi(k - 1, :) .* xip(k - 1, :);
        X(k, :) = (f(k, :) - xi(k - 1, :) .* X(k - 1, :)) ./ w;
    end

    for k = K - 1:-1:1
        X(k, :) = X(k, :) - xip(k, :) .* X(k + 1, :);
    end
