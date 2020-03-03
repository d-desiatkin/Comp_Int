pkg load optim
pkg load geometry

%ext_f = [rand([2 1])*2.5; -10];
%ext_t = rand([3 1])*2.5;
ext_f = [0; 0; 10];
ext_t = [0; 0; 0];

f_init = [0.07; 0.03; 0.08; 0.01; 0.001; 0.001; 0.001; 0.001; 0.001];

% Robot legs radius vectors %
rk1 = [0; 2/3; 0];
rk2 = [-1/sqrt(3); -1/3; 0];
rk3 = [1/sqrt(3); -1/3; 0];

% Equality constraints %
Aeq = [eye(3) eye(3) eye(3);
       skew(rk1) skew(rk2) skew(rk3)];
Beq = [ext_f;
       ext_t];
       
% Inequality constraints %
% Lets choose the cone with cirle radius 1 and height 1 %
% Cone vectors in leg reference frame %
rl1 = [0;1;1];
rl2 = [-1;0;1];
rl3 = [0;-1;1];
rl4 = [1;0;1];

cv_1 = normalizeVector(cross(rl2, rl1));
cv_2 = normalizeVector(cross(rl3, rl2));
cv_3 = normalizeVector(cross(rl4, rl3));
cv_4 = normalizeVector(cross(rl1, rl4));

cone1 = [cv_1 cv_2 cv_3 cv_4]';

A = [coneGen([0;-1;1], 1, [], 0, 8) zeros(8,6);
         zeros(8,3) coneGen([0;0;1], 1, [], 0, 8) zeros(8,3);
         zeros(8,6) coneGen([0;0;1], 1, [], 0, 8)];
         
B = zeros(24,1);


H = eye(9);
F = zeros(9,1); 

[forces, ~, flag] = quadprog(H, F, A, B, Aeq, Beq, [], []);

forces, flag

A*forces

pkg unload geometry
pkg unload optim
