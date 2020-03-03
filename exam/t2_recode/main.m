ext_f = [0; 0; -10];
ext_t = [0; 0; 0];

rk1 = [0; 2/3; 0];
rk2 = [-1/sqrt(3); -1/3; 0];
rk3 = [1/sqrt(3); -1/3; 0];

% Equality constraints %
F = [eye(3) eye(3) eye(3);
       skew(rk1) skew(rk2) skew(rk3)];
g = [-ext_f;
     ext_t];

th = pi/4;
fT = ones(9,1)';
mu = 0.5;


A = eye(9);
d = [0; 0; 1];
%A = [skew(d) skew(d) skew(d);
%        d' d' d'];
B = -ext_f'*d*mu;

cvx_begin
variable forces(9,1)
minimize fT*forces
    subject to
       F*forces == g;
       A*forces <= B;

cvx_end
