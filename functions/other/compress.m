function y = compress(x, factor)
    cut = [1:factor:numel(x), numel(x)];
    y = zeros(numel(cut)-1, 1);
    for idx = 1:numel(cut)-1
        y(idx) = mean(x(cut(idx):cut(idx+1)));
    end
end

