function [M_mtrx, n_vctr] = dynamics_mat(M,m,L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

syms x x_d x_2d thta thta_d thta_2d         real

g0 = 9.81;

q = [x, thta]';
q_d = [x_d, thta_d]';
q_2d = [x_2d, thta_2d]';

x_p = x + L*sin(thta);
y_p = L*cos(thta);

x_p_d = jacobian(x_p,q)*q_d;
y_p_d = jacobian(y_p,q)*q_d;

K_cart = 0.5*M*x_d^2; 
K_pend = simplify(0.5*m*(x_p_d^2 + y_p_d^2));
K = K_cart + K_pend;

P_pend = m*g0*y_p;

Lgr = K - P_pend;

t1 = jacobian(Lgr,q); % dL/dq
t2 = jacobian(Lgr,q_d);% dL/dq_d
t3 = jacobian(t2, [q;q_d])*[q_d;q_2d]; % d/dt (dL/dq_d)
t4 = simplify(t3-t1'); % t4 = M(q)*q_2d + n(q,q_d)
clear t1 t2 t3 

M_mtrx = jacobian(t4,q_2d);
n_vctr = simplify(t4 - M_mtrx*q_2d);
B = [1 0]';

matlabFunction(M_mtrx,'File', 'M_mtrx','Vars', {q});
matlabFunction(n_vctr,'File', 'n_vctr','Vars', {q,q_d});

end

