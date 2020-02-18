% % g = admDiffFor(@myFun, 1, x)
% % 
% % H = admHessFor(@myFun, 1, x)
% % 
% % p = -H\g'

x = linspace(-10,10,100);
y = linspace(-10,10,100);
[x,y] = meshgrid(x,y);
contour3(x,y,x.^4+y.^2,100)
%axis('equal');
view(0,90)
hold on


%% redefine function for optimization.
% in the form that it accepts vector input
f = @(v)(elliplticParaboloid(v));
df = @(v)(admDiffFor(@elliplticParaboloid, 1, v)');


x0 = [10;10];
tol = 1e-8;
err = 1;

for k = 1:1000
     %g = admDiffFor(@elliplticParaboloid, 1, x0);
     %H = admHessFor(@elliplticParaboloid, 1, x0);
     %dir = -H\g';
        
    dir = -df(x0)
        
    alpha = lineSearchWolfe(f,df,x0,dir)
    %alpha = mb_nocLineSearch(f,df,x0,dir,df(x0)'*dir,f(x0))
    
    x1 = x0 + alpha*dir;
    
    plot([x0(1),x1(1)],[x0(2),x1(2)],'--o' );
    pause(0.01)
    
    err = norm(x1-x0);
    if err < tol
        break
    end    
    
    x0 = x1;
    %alpha = alphaStart;
    
end