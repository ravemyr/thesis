    function t = H(mesh,points,q,xi,eta,Q)
N = max(size(points));
d = min(size(points));
if(d==2)
    M = sqrt(max(size(mesh)));
    f = zeros(M,M);
    for i = 1:N
        for ii = 1:M
            for iii = 1:M
                %for j = -Q/2:Q/2
                    f(ii,iii) = f(ii,iii) + q(i)*((-2*xi^2)/(pi*eta))^(3/2) ...
                    *exp((-2*xi^2*norm(mesh(8*(ii-1)+iii,:)-points(i,:))^2)/eta);
                %end
            end
        end
    end
elseif(d==1)
    M = max(size(mesh));
    f = zeros(M,1);
    
    for i = 1:N
       [~,pos] = min(abs(mesh(:,:)-points(i,:)));
       for ii = -(Q-1)/2:(Q-1)/2
           if(pos+ii==M+1)
               pos = pos-M;
           elseif(pos+ii<1)
               pos = pos + M;
           end
            f(pos+ii) = f(pos+ii) + q(i)*((2*xi^2)/(pi*eta))^(3/2) ...
            *exp(-2*xi^2*(norm(mesh(pos+ii,:)-points(i,:))^2)/eta);
       end
    end
end
t = f;