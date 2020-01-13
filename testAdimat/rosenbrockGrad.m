function r = rosenbrockGrad(x)
r = [-400*(x(2) - x(1).^2) - 2*(1-x(1));
     200*(x(2) - x(1).^2)];
 