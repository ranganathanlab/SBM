function [W]=Wj(J,h)
% convert from Coupling and Field matrices parameters to one long vector of parameters
% note that since we work in the Ising Gauge, some of the variables are still redundant
    q=size(J,1);
    n=size(J,3);
    W=zeros(q*n+q*q*n*(n-1)/2,1);
    count=1;
    for i=1:n
        for j=i+1:n
            for a=1:q
                for b=1:q
                    W(count)=J(a,b,i,j);
                    count=count+1;
                end
            end
        end
    end
    for i=1:n
        for a=1:q
            W(count)=h(a,i);         
            count=count+1;
        end
    end
end