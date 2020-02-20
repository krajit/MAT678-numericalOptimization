function algo_4p1
clc
x = linspace(-10,10,100);
y = linspace(-10,10,100);
[x,y] = meshgrid(x,y);


subplot(1,2,1)
contour3(x,y,ellipticContour(x,y),100)
%axis('equal');
view(0,90)
hold on

%% redefine function for optimization.
% in the form that it accepts vector input
f = @(v)(elliptic(v));
df = @(v)(admDiffFor(@elliptic, 1, v)');
Hf = @(v)(admHessFor(@elliptic, 1, v)');


x0 = [10;10];
h0 = [0;0];
tol = 1e-10;
err = 1;
tau_k = 1;
Delta_k = 1;
Delta_bar = 2;
eta = 1/6;

for k = 1:1000
  
    % find gradient and hessian
    gk = df(x0);
    Bk = Hf(x0);
    
    
    % create mk    
    mk = @(h)(f(h0) + gk'*h + 0.5*h'*Bk*h);
    
    % get pk
    %pk = CauchyPk();
    pk = dogLegPk();
    
    
    % evaluate rho_k
    rho_k = (f(x0) - f(x0 + pk))/(mk(h0) - mk(pk));
    
    if rho_k < 0.25
        Delta_k = 0.25*norm(pk);
    else
        if rho_k > 0.75 && abs(norm(pk) - Delta_k) < 1e-8
            Delta_k = min(2*Delta_k, Delta_bar);
        end
    end
    
    if rho_k > eta
        x1 = x0 + pk;
    else
        x1 = x0;
    end
    
    subplot(1,2,1)
    plot([x0(1),x1(1)],[x0(2),x1(2)],'--o' );
    pause(0.01)
    
    subplot(1,2,2)
    plot([k-1 k],[f(x0) f(x1)])
    hold on
    
%     err = norm(x1-x0);
%     if err < tol
%         break
%     end    
    k
    x0 = x1
end

%% get Cauchy pk
    function pk = CauchyPk()
            % get tau_k
        if (gk'*Bk*gk <= 0)
            tau_k = 1;
        else
            tau_k = min(norm(gk)^3/(Delta_k*(gk'*Bk*gk)),1);
        end
    
        % get p_k
        pk = -tau_k*(Delta_k/norm(gk))*gk;
    end


%% get dog-leg pk
    function pk = dogLegPk()
        pu = -((gk'*gk)/(gk'*Bk*gk))*gk;
        pb = -Bk\gk;
        
        if (norm(pu)>= Delta_k)
            pk = Delta_k*pu/norm(pu);
        else
            a = norm(pb-pu)^2;
            b = 2*sum(pu.*(pb-pu));
            c = norm(pu)^2 - Delta_k^2;
            
            s1 = (-b+sqrt(b^2 - 4*a*c))/(2*a);
            s2 = (-b-sqrt(b^2 - 4*a*c))/(2*a);
            t1 = s1+1;
            t2 = s2+1;
            
            t = min(2,max(1,t1));            
            pk = pu + (t-1)*(pb-pu);           
            
        end           
      
    end

end
