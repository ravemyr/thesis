function outopt =  param_select_stokes(tol,inopt)
    %Control that at least L and xi are given, compute x and f if not
    %given, select N =10 if not given
    checkInopt(inopt);
    outopt.xi = inopt.xi;
    F = sum(norm((inopt.f).^2));
    B = sqrt(F)*inopt.xi;
    est_in = 8*F/(3*(pi*inopt.box(1)*tol)^2);
    M = ceil(2*sqrt((inopt.xi*inopt.box(1))^2/(2*pi^2)*lambertw(est_in)));
   
%     if(lim(xi*L/M)>val)
%        %TODO Modify selection of parameters if threshold not met 
%     end
    if ~isfield(inopt,'window'), inopt.window = 'gaussian'; outopt.window=inopt.window; end
    if(strcmp(inopt.window,'gaussian'))
        %Compute P from error estimate here
        c = sqrt(0.91);
        outopt.P = ceil(-log(tol/(2*B))/(2*pi*c));
        outopt.w = outopt.h*outopt.P/2;
        outopt.m = 0.95*sqrt(pi*outopt.P);  
        outopt.eta = (2*outopt.w*outopt.xi/outopt.m)^2;
        outopt.c = 2*outopt.xi^2/outopt.eta;
    end
    if strcmp(inopt.window,'kaiser_exact') ...
        || strcmp(inopt.window,'kaiser_poly')
      outopt.betaP = 2.5;   
      fx = @(x) 0.30*x^2 + 1.73 * x - 1.62;
      disp(M)
      if(fx(M/(inopt.xi*inopt.box(1)))>ceil(-log(tol/(10*B))/2.5))
            outopt.P = ceil(-log(tol/(10*B))/2.5);
      else
	if(mod(M,2)==0) M = M+12; else M = M+13; end
        outopt.P = ceil(-log(tol/(10*B))/2.5)+16;
      end
      outopt.beta = outopt.betaP*outopt.P;
      outopt.kaiser_scaling = 1/besseli(0,outopt.beta);
    end
    if strcmp(inopt.window,'kaiser_poly')
      outopt.polynomial_degree = min(outopt.P/2 +2,9);
    end
    if mod(outopt.P,2)==0
      outopt.p_half = outopt.P/2;
    else
      outopt.p_half = (outopt.P-1)/2;
    end
    outopt.M = [M,M,M];
    outopt.box = inopt.box;
end
    %Verify that user data exists, generate if not
function checkInopt(opt)
    assert(isfield(opt, 'box'), 'cell size box must be given in opt struct');
    %Assume xi given
    assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
    if(~isfield(opt,'x'))
        if(~isfield(opt,'N')) opt.N = 10; end
        [opt.x, opt.f] = SE_charged_system(opt.N,opt.box,'vector'); 
    end
    % Half-support of window
end
