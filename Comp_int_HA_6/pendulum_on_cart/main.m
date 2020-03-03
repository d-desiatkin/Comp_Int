clc
clear all
close all

% System parameters
M = 1;
m = 0.3;
L = 0.5;
% System dynamics
[M_mtrx, n_vctr] = dynamics_mat(M,m,L);

% Number of points
n = 30;

% boundary constraints
d_max = 2; % max cart shift in both directions
t_max = 2; % end time

x_0 = [0; pi; 0; 0; 0]; % starting conditions
x_f = [1; 0; 0; 0; 0]; % end conditions

% parameters creation for trapezoidal collocations
time = linspace(0, t_max, n);
x = [sin(linspace(x_0(1), x_f(1)+2*pi, n));
     linspace(x_0(2), x_f(2), n);
     linspace(x_0(3), x_f(3), n);
     linspace(x_0(4), x_f(4), n);
     linspace(x_0(5), x_f(5), n)];
x = reshape(x, [], 1);


% % Filling equality constraints trapezoidal collocations
% Aeq = zeros(10, 5*n);
% Aeq(1:5, 1:5) = eye(5);
% Aeq(6:10, 5*n-4: 5*n) = eye(5);
% beq = [x_0;x_f];
% 
% options = optimoptions("fmincon", "Display", "iter", "MaxFunctionEvaluations", 40000);
% my_con = @(x)mycon_trap(time, x);
% my_cost = @(x)costFun_trap(time, x);
% x = fmincon(my_cost, x, [],[],Aeq,beq,[],[], my_con, options);
% tmp = reshape(x, 5, n);
% min_dec_var = min(tmp, [], 2);
% max_dec_var = max(tmp, [], 2);
% sprintf("min / max x: %f, %f",min_dec_var(1), max_dec_var(1))
% sprintf("min / max theta: %f, %f",min_dec_var(2), max_dec_var(2))
% sprintf("min / max dx: %f, %f",min_dec_var(3), max_dec_var(3))
% sprintf("min / max dtheta: %f, %f",min_dec_var(4), max_dec_var(4))
% sprintf("min / max control input: %f, %f",min_dec_var(5), max_dec_var(5))
% visualize_cart_pendulum(x);


% parameters creation for Hermite - Simpson collocations
time = linspace(0, t_max, 2*n+1);
x = [sin(linspace(x_0(1), x_f(1)+2*pi, 2*n+1));
     linspace(x_0(2), x_f(2), 2*n+1);
     linspace(x_0(3), x_f(3), 2*n+1);
     linspace(x_0(4), x_f(4), 2*n+1);
     linspace(x_0(5), x_f(5), 2*n+1)];
x = reshape(x, [], 1);


% Filling equality constraints Hermite - Simpson collocations
Aeq = zeros(10, 5*(2*n+1));
Aeq(1:5, 1:5) = eye(5);
Aeq(6:10, 5*(2*n+1)-4: 5*(2*n+1)) = eye(5);
beq = [x_0;x_f];

options = optimoptions("fmincon", "Display", "iter", "Algorithm", "sqp");
my_con = @(x)mycon_HS(time, x);
my_cost = @(x)costFun_HS(time, x);
[x, fval] = fmincon(my_cost, x, [],[],Aeq,beq,[],[], my_con, options);
tmp = reshape(x, 5, 2*n+1);
min_dec_var = min(tmp, [], 2);
max_dec_var = max(tmp, [], 2);
output = [
sprintf("min / max x: %f, %f",min_dec_var(1), max_dec_var(1));
sprintf("min / max theta: %f, %f",min_dec_var(2), max_dec_var(2));
sprintf("min / max dx: %f, %f",min_dec_var(3), max_dec_var(3));
sprintf("min / max dtheta: %f, %f",min_dec_var(4), max_dec_var(4));
sprintf("min / max control input: %f, %f",min_dec_var(5), max_dec_var(5));
sprintf("loss function: %f", fval)]
figure(1);
hold on;
visualize_cart_pendulum(x);
hold off;
figure(2);
hold on;
plot(tmp(5, :));
