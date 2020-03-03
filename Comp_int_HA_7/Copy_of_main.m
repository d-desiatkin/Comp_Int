clc
clear all
close all

% SOLVING ODE

optns = odeset('RelTol',1e-6,'AbsTol',1e-6,'NormControl','on');
tspan = [0, 1];
x_0 = [0; 0; 0.015; 0];
[t,x] = ode45( @(t,x)dynamics(t,x),tspan,x_0,optns);


% x3max = ;
x1max = 0.02;
tau = 0.4;
% R = diag([inf; inf; 1/x3max^2; inf]);
x1const = zeros(length(t),1);
SWN = zeros(length(t),1);
Energy = zeros(length(t),1);
for i=1:length(t)
    G = diag([1/(x1max*exp(-t(i,1)/tau))^2; 0; 0; 0]);
    x1const(i,1) = x1max*exp(-t(i,1)/tau);
    SWN(i,1) = x(i,:)*G*(x(i,:)');
end

f1 = figure();
plot(t,SWN);

f2 = figure();
plot(t, x(:,1));
hold on;
plot(t,x(:,3), 'g--');
plot(t, x1const,'r-.', t, -x1const, 'r-.');
legend('x1', 'x3', 'x1_{constr.}')
hold off;


Ks = 50000;
Ku = 200000;
Ms = 320;
Mu = 40;
Bs = 2000;

A = [0 1 0 -1
    -Ks/Ms -Bs/Ms 0 Bs/Ms
     0 0 0 1
     Ks/Mu Bs/Mu -Ku/Mu -Bs/Mu];

x1max = 0.02;
tau = 0.4;
x3max = 0.0151;
eps = 1e-6;
N = length(t);

R = diag([1/eps; 1/eps; 1/x3max^2; 1/eps]);
G = diag([1/(x1max*exp(-t(1,1)/tau))^2; eps; eps; eps]);

cvx_begin sdp
cvx_precision best
    variable P(4,4) symmetric
    variable P_dot(4,4)

    subject to
        P - R < 0;
        P - G > 0;
        P_dot + A'*P + P*A < 0;
cvx_end

% Сheck first point conditions.
[~,p1] = chol(-(P_dot + A'*P + P*A));
[~,p2] = chol((P - G));
[~,p3] = chol(-(P - R));
p = [p1;p2;p3];
if all(p == 0)
    fprintf("Conditions was satisfied for first point\n\n")
else
    fprintf("Conditions wasn't satisfied for first point\n\n")
    
end

P_first = P;
P_dot_prev = P_dot;
P_cum = zeros(4,4);
p1 = zeros(N-1,1);
p2 = zeros(N-1,1);

for i=2:N
    
    G = diag([1/(x1max*exp(-t(i,1)/tau))^2; eps; eps; eps]);
    
    
    cvx_begin sdp
        variable P_dot(4,4)
        variable P(4,4) symmetric
        subject to
            P == P_first + P_cum + P_dot_prev * (t(i,1)-t(i-1,1));
            P - G > 0;
            P_dot + A'*P + P*A < 0;
    cvx_end
    [~, p1(i-1,1)] = chol(-(P_dot + A'*P + P*A));
    [~, p2(i-1,1)] = chol((P - G));
    P_cum = P_cum + P_dot_prev * (t(i,1)-t(i-1,1));
    P_dot_prev = P_dot;

end

p = [p1; p2];

if all(p == 0)
    fprintf("Conditions was satisfied for consequative points\n\n")
else
    fprintf("Conditions wasn't satisfied for consequative points\n\n")
end


