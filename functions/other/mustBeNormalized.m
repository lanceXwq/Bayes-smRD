function mustBeNormalized(M)

    if any(sum(M, 2) - 1 > eps)
        error('This matrix is not normalized.')
    end

end
