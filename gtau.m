function res = gtau(data,t)
    N = length(data);
    res = zeros(1,N);
    p = 1;
    for x = data
        for q = -1000:1000
            res(p) = res(p) + exp((-1/(4*t))*(x-2*pi*q)^2);
        end
        p = p+1;
    end
end