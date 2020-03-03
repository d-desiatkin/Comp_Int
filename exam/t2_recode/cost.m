function [cost] = cost(x)
%COST Summary of this function goes here
%   Detailed explanation goes here
f = [1, 1, 1, 1, 1, 1, 1, 1, 1];
cost = 0;
tmp = reshape(x, 9,1);

ext_f = [0; 0; -10];
ext_t = [0; 0; 0];

rk1 = [0; 2/3; 0];
rk2 = [-1/sqrt(3); -1/3; 0];
rk3 = [1/sqrt(3); -1/3; 0];


% Equality constraints %
F = [eye(3) eye(3) eye(3);
       skew(rk1) skew(rk2) skew(rk3)];
g = [-ext_f;
     -ext_t];

th = pi/4;
fT = ones(9,1)';
mu = 0.5;


A = eye(9);
d = [0; 0; 1];
%A = [skew(d) skew(d) skew(d);
%        d' d' d'];
B = -ext_f'*d*mu;
eps = 10;

x1 = (A*tmp - B);
y11 = norm(normpdf(x1(1:3,:), 0, eps));
y12 = norm(normpdf(x1(4:6,:), 0, eps));
y13 = norm(normpdf(x1(7:9,:), 0, eps));

cost = cost - (y11 + y12 + y13);


x2 = F*tmp - g;
y21 = norm(normpdf(x2(1:3,:), 0, eps));
y22 = norm(normpdf(x2(4:6,:), 0, eps));


cost = cost - (y21+y22)*10;


end

