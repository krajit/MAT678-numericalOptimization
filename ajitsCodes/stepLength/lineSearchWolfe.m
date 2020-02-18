function aS = lineSearchWolfe(f,df,xk,pk)

% f    - objective function handle
% df  - gradient function handle
% xk - Current value
% pk - search direction

% define restriction of f & f' on the line xk + a*pk
phi = @(a)(f(xk + a*pk));
phip = @(a)(df(xk+a*pk)'*pk);

c1 = 1e-4;
c2 = 0.8;

aMax = 100;
aim1 = 0;
ai = 1;
iter = 0;


while 1
    if (phi(ai) > phi(0)+c1*ai*phip(0)) || (phi(ai) >= phi(aim1))
        aS = zoom(aim1,ai);
        return
    end
    
    if (abs(phip(ai)) <= -c2*phip(0))
        aS = ai;
        return
    end
    
    if phip(ai) >= 0
        aS = zoom(ai,aMax);
        return
    end
    
    ami1 = ai;
    ai = min(aMax,ai*3); % idea copied from MB;    
    
    iter = iter + 1;
end

    function aj = zoom(aLo, aHi)
        nMax = 100;  % max iteration in zooming        
            for k = 1:nMax
                aj = (aLo + aHi)/2; % using bisection now. 
                                                % TODO. Cubic, quadratic
                if (phi(aj) > phi(0) + c1*aj*phip(0)) || (phi(aj) >= phi(aLo))
                    aHi = aj;
                else
                    if (abs(phip(aj)) <= -c2*phip(0))
                        return
                    end

                    if phip(aj)*(aHi - aLo) >= 0
                        aHi = aLo;
                         aLo = aj;
                    end
                   

                end
            end  
        end
end