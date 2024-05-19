function X =bits_to_4pam(b1,b2)
    N=100;
    for i=1:N/2
           if (b1(i)==0 & b2(i)==0)
               X(i)=3;
           elseif (b1(i)==0 & b2(i)==1)
               X(i)=1;
           elseif (b1(i)==1 & b2(i)==1)
               X(i)=-1;
           else 
               X(i)=-3;
           end
       end
end