function outopt =  param_select_stokes(tol,inopt)%,ErrorType
    %Control that at least L and xi are given, compute x and f if not
    %given, select N =10 if not given
%     tol = varargin{1};
%     inopt = varargin{2};
%     ErrorType = varargin{3};
%     checkInopt(inopt);
    outopt.xi = inopt.xi;
    F = sum(norm((inopt.f).^2));
    B = sqrt(F)*inopt.xi;
%     if(strcmp(ErrorType,'Relative'))
%        tol = tol/B;
%     end
    est_in = 8*F*inopt.xi^2/(3*(pi*tol)^2);
    M = ceil(2*sqrt((inopt.xi*inopt.box(1))^2/(2*pi^2)*lambertw(est_in)));
    if ~isfield(inopt,'window'), inopt.window = 'gaussian'; outopt.window=inopt.window; end
%     if(strcmp(inopt.window,'gaussian'))
%         %Compute P from error estimate here
%         c = sqrt(0.91);
%         outopt.P = ceil(-2*log(tol)/(pi*c)) + 11;
%         if(inopt.xi*inopt.box(1)>30)
%             if(mod(ceil(M),2)==0) M = ceil(M); else M = ceil(M); end
%             M = M+ceil(0.75*inopt.xi*inopt.box(1));
%             %M = M + ceil(-log10(tol)/2);
%         else	
%         M = M + ceil(0.5*inopt.xi*inopt.box(1));%+max(ceil(-log10(tol)/2-3),0);
%         end end
    if (strcmp(inopt.window,'kaiser_exact') || strcmp(inopt.window,'kaiser_poly'))
	outopt.P = ceil(-log(tol/10)/2.5)+2;
	M = ceil(1.18*M);
    disp(M)  	
    end
    if strcmp(inopt.window,'kaiser_poly')
      outopt.polynomial_degree = min(ceil(outopt.P/2 +2),9);
      if(mod(outopt.P,2)~=0&&outopt.P<16)
        outopt.P = outopt.P +1;
      else
        outopt.P = 16;
	disp('Kaiser poly supports upto P=16 in this version')
      end	
    end



   
    outopt.window = inopt.window;
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
        %[opt.x, opt.f] = SE_charged_system(opt.N,opt.box,'vector'); 
    end
    % Half-support of window
end
