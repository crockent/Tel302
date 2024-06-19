%12)
function est_X = detect_4_PAM(Y, A)
    levels = [-3 -1 1 3] * A;
    est_X = zeros(size(Y));
    for i = 1:length(Y)
        [~, idx] = min(abs(levels - Y(i)));
        est_X(i) = levels(idx);
    end
end