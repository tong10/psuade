function Y = ish(X1,X2,X3)

X1a = X1 * 2 * pi - pi;
X2a = X2 * 2 * pi - pi;
X3a = X3 * 2 * pi - pi;
sin(X1a)
7.0*sin(X2a)*sin(X2a)
0.1*X3a*X3a*X3a*X3a*sin(X1a)

Y = sin(X1a)+7.0*sin(X2a)+0.1*X3a*X3a*X3a*X3a*sin(X1a);

return

