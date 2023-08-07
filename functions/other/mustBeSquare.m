function mustBeSquare(M)
    s = size(M);

    if length(s) ~= 2 || s(1) ~= s(2)
        error('This is not a square matrix.')
    end

end
