    function t = H(mesh,points,q,xi,eta,Q)
N = max(size(points));
M = sqrt(max(size(mesh)));
f = zeros(M,M);
for i = 1:N
    for ii = 1:M
        for iii = 1:M
            %for j = -Q/2:Q/2
                f(ii,iii) = f(ii,iii) + q(i)*((2*xi^2)/(pi*eta))^(3/2) ...
                *exp((-2*xi^2*norm(mesh(8*(ii-1)+iii,:)-points(i,:))^2)/eta);
            %end
        end
    end
end
t = f;