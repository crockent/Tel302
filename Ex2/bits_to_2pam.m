function [b] =bits_to_2pam(b)
    for i=1:length(b)
        if b(i)==0
            b(i)=1;
        elseif b(i)==1
            b(i)=-1;
        end
    end
end