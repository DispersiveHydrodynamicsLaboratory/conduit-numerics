function [ Y ] = connect_fit( x , parameters )

if nargin < 2
    parameters = [ 1 2 1 2 1 0];
end

x1    = parameters(1);     x2    = parameters(2);
y1    = parameters(3);     y2    = parameters(4);
dydx1 = parameters(5);     dydx2 = parameters(6);

a = x(1); b = x(2); c = x(3); d = x(4);

X1 = a*log(cosh(b*(x1-c))) + d - y1;
X2 = a*log(cosh(b*(x2-c))) + d - y2;
X3 = a*b*tanh(b*(x1-c)) - dydx1;
X4 = a*b*tanh(b*(x2-c)) - dydx2;

Y = [X1;X2;X3;X4];
end

