function outopt =  param_select_stokes(tol,inopt)
    %Control that at least L and xi are given, compute x and f if not
    %given, select N =10 if not given
    checkInopt(inopt);
    outopt.xi = inopt.xi;
    F = sqrt(sum(norm((inopt.f).^2)));
    B = F*inopt.xi;
    kf_inside = lambertw(4*(inopt.box(1)*F)^(2/3)*inopt.xi^2/(3*pi^(10/3)*tol^(4/3)));
    M = ceil(sqrt(3)*inopt.box(1)*inopt.xi*sqrt(kf_inside)/pi);
   outopt.h = M./inopt.box(1); 
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
      fx = @(x) 0.20*x^2 + 2.3 * x - 2.2527;
      if(fx(M/(inopt.xi*inopt.L))>ceil(-log(tol/(10*B))/2.5))
      outopt.P = ceil(-log(tol/(10*B))/2.5)+5;
      else 
          outopt.M = outpot.M + 8;
          outopt.P = ceil(-log(tol/(10*B))/2.5)+5;
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
