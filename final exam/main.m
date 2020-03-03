clc
clear all
close all

l = [1, 0.7];
lcm = l/2;
mass = [10, 7];
g = 0;

I_1 = 1 / 12 * mass(1) * lcm(1) ^ 2;
I_2 = 1 / 12 * mass(2) * lcm(2) ^ 2;

syms q1 q2 x y q1d q2d xd yd q12d q22d x2d y2d   real

q = [q1 q2]';
q_d = [q1d q2d]';
q_2d = [q12d q22d]';

T = Rz(q1)*Tx(l(1))*Rz(q2)*Tx(l(2));
T1 = Rz(q1)*Tx(lcm(1));
T2 = Rz(q1)*Tx(l(1))*Rz(q2)*Tx(lcm(2));

vel1 = jacobian(T1(1:2,4), q)*q_d;
vel2 = jacobian(T2(1:2,4), q)*q_d;

K = 1/2*(vel1'*mass(1)*vel1 + vel2'*mass(2)*vel2);

P = g*(mass(1)*T1(2,4) + mass(2)*T2(2,4));

L = K - P;

t1 = jacobian(L,q); % dL/dq
t2 = jacobian(L,q_d);% dL/dq_d
t3 = jacobian(t2, [q;q_d])*[q_d;q_2d]; % d/dt (dL/dq_d)
t4 = simplify(t3-t1'); % t4 = M(q)*q_2d + n(q,q_d)
clear t1 t2 t3

M_mtrx = jacobian(t4,q_2d);
n_vctr = simplify(t4 - M_mtrx*q_2d);
clear t4

matlabFunction(M_mtrx,'File', 'M_mtrx','Vars', {q});
matlabFunction(n_vctr,'File', 'n_vctr','Vars', {q,q_d});

% SOLVING ODE
tspan = [0, 2];
x_0 = [60,-30,0,0]'*pi/180;

cost([1,1,1,1],tspan,x_0);


lb = ones(1,4)*(-10);
ub = ones(1,4)*(10);
params = [sprintf("Initialize particle swarm");...
          sprintf("Lower bounds: %f %f %f %f ", lb);...
          sprintf("Upper bounds: %f %f %f %f ", ub)];
disp(params);
ans = particleswarm(@(x)cost(x, tspan, x_0), 4, lb, ub);
% bbb = reshape(ans(1:6*aaa), 6, aaa);
K = reshape(ans, 2, 2);

% visualization(bbb, l);







