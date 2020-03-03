pkg load optim
pkg load geometry


% Parameters
n = 100
T = 1
dt = T/n
mass = 1;

% Understanding the matix sizes
x0 = zeros(n*2, 1);
Aeq = zeros(2+2+n, n*2);
Beq = zeros(2+2+n, 1);
A = zeros(n, 2*n);
B = zeros(n, 1);
F = zeros(n*2,1);
H = 0.5 * [zeros(n, 2*n); zeros(n,n), eye(n)/mass^2];

% Filling the equality constraints
Aeq(1, 1:2*n) = [1, zeros(1, 2*n-1)];
Aeq(2, 1:2*n) = [zeros(1, n-1), 1, zeros(1, n)];
Aeq(3, 1:2*n) = [-1/dt, 1/dt, zeros(1, n-2), zeros(1 ,n)];
Aeq(4, 1:2*n) = [zeros(1, n-2), -1/dt, 1/dt, zeros(1, n)];
%disp(Aeq(1:4,:));
Aeq(5, 1:3) = [1 -2 1];
Aeq(5, n+1) = -dt^2/mass;
Aeq(6:n+3,1:2*n) = [eye(n-2) zeros(n-2, n+2)] -...
               2*[zeros(n-2,1) eye(n-2) zeros(n-2, n+1)] +...
               [zeros(n-2,2) eye(n-2) zeros(n-2, n)] -...
               dt^2/mass*[zeros(n-2,n+1) eye(n-2) zeros(n-2, 1)];
Aeq(n+4, n-2:n) = [1 -2 1];
Aeq(n+4, 2*n) = -dt^2/mass;

Beq(3,1) = 1;
Beq(4,1) = -1;

% Filling the inequality constraints
A = [eye(n), zeros(n,n)];
B = B + 1/9;
x0 = [0; ones(n-2,1)*1/18; 0; -ones(n,1)];


% Solver
options = optimset('MaxIter', 2000);
[x, J, flag] = quadprog(H, F, A, B, Aeq, Beq, [], [], x0, options);

%A * coord
flag
J
A*x
t = 0:T/(n-1):1;
plot(t, x(n+1:2*n,1));
xlabel ("time (s)");
ylabel ("force (N)");
title ("Applied force");

pkg unload geometry
pkg unload optim
