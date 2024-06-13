%%2%%
function X =bits_to_4_PAM(bits,A)
    X = zeros(1, length(bits) / 2);
    grey_code = [-3*A,-A,A,3*A];
    k=1;
    l=2;
    for i= k: length(bits)/2    
        if(bits(k)==0 & bits(l)==0)
            X(i)=grey_code(1);
        elseif (bits(k)==0 & bits(l)==1)
            X(i)=grey_code(2);
        elseif (bits(k)==1 & bits(l)==1)
            X(i)=grey_code(3);
        else 
            X(i)=grey_code(4);
         end
        k = k+2;
        l = l+2
    end
end