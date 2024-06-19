%14)
function est_bit = PAM_4_to_bits(X, A,b)
    levels = [-3 -1 1 3] * A;
    est_bit = zeros(1,length(b));

 for i = 1:length(X)
        if X(i) == levels(1)
            est_bit(2*i-1:2*i) = [0 0];
        elseif X(i) == levels(2)
            est_bit(2*i-1:2*i) = [0 1];
        elseif X(i) == levels(3)
            est_bit(2*i-1:2*i) = [1 1];
        else % X(i) == levels(4)
            est_bit(2*i-1:2*i) = [1 0];
        end
    end
end
