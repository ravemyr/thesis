function res = gtau(data,t,sp)
    N = length(data);
    res = zeros(1,N);
    p = 1;
    for x = data
        for q = -sp:sp
            res(p) = res(p) + exp((-1/(4*t))*(x-2*pi*q)^2);
        end
        p = p+1;
    end
end