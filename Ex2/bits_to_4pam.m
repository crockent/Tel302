function X =bits_to_4pam(b1,b2)
    
    for i=1:length(b1)
       for j=1:length(b2)
           if (b1(i)==0 & b2(j)==0)
               X(i)=3;
           elseif (b1(i)==0 & b2(j)==1)
               X(i)=1;
           elseif (b1(i)==1 & b2(j)==1)
               X(i)=-1;
           elseif (b1(i)==1 & b2(j)==0)
               X(i)=-3;
           end
       end
end