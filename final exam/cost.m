function [cost] = cost(K, tspan, x_0)
%COST Summary of this function goes here
%   Detailed explanation goes here

K = reshape(K, 2, 2);

optns = odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on');
[t, x] = ode45( @(t,x)dynamics(t, x, K, x_0),tspan, x_0, optns);
[~, traj, u] = dynamics(t, x, K, x_0);

[n,d] = size(x);

cost = norm(K, 'fro')*100;



for i=1:n
    cost = cost + (x(i,1:2) - traj(i,1:2))*eye(2)*(x(i,1:2)- traj(i,1:2))';
    cost = cost + u(i,1:2)*eye(2)*u(i,1:2)';
end

% Equality constraints
for i = 1:n-1
    f_1 = dynamics(t(i,:)', x(i,:)', K, x_0);
    f_2 = dynamics(t(i+1,:)', x(i+1,:)', K, x_0);
    rew_tmp = (t(i+1)-t(i))/2 * (f_2 + f_1) + x(i , 1:4)' - x(i+1, 1:4)';
    reward = -(100 - rew_tmp'*eye(4)*rew_tmp);
    cost = cost + reward;
end

% Inequality constraints
A1 = diag([1,1,0,0]);
B1 = repmat([1;1;1;1], n, 1)*pi/180*90;
B2 = repmat([1;1;1;1], n, 1)*pi/180*0;

if all(reshape(A1*x', n*4, 1) <= B1)
    cost = cost + sum(-log(B1 - reshape(A1*x', n*4, 1)));
    disp(K);
    disp(cost); 
    disp(norm(K, 'fro'));
else
    cost = Inf;
end

if all(-reshape(A1*x', n*4, 1) <= B2)
    cost = cost + sum(-log(reshape(A1*x', n*4, 1)-B2) );
    disp(K);
    disp(cost); 
    disp(norm(K, 'fro'));
else
    cost = Inf;
end



end