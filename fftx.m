function y = fftx(x)
    n = max(size(x));
    omega = exp(-2*1i*pi/n);
    if(mod(n,2)==0)
        k = (0:(n/2-1))';
        w = omega.^k;
        %Special case management
        if(x(1:2:n-1)==x(2:2:n))
            u = fftx(x(1:2:n-1));
            v = w.*u;
        else
            %Recursive step
            u = fftx(x(1:2:n-1));
            v = w.*fftx(x(2:2:n));
        end
        y = [u+v;u-v];
    else
        %Base step
        j = 0:n-1; k = j';
        F = omega.^(k*j);
        y = F*x';
    end
end