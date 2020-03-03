function [c,ceq] = mycon_HS(t, x)
%MYCON Summary of this function goes here
%   Detailed explanation goes here
c = [];  % Compute nonlinear inequalities at x.
ceq = [];   % Compute nonlinear equalities at x.

l = length(x);
x = reshape(x, 5, l/5);
[~, n] =  size(x);

ceq = zeros(4, n-2);

for i=1:2:n-2
    h_k = t(i+2)-t(i);
    x_1 = x(:, i);
    x_2 = x(:, i+2);
    f_1 = dynamics(t(i), x_1);
    f_2 = dynamics(t(i+2), x_2);
    x_m = x(:, i+1);
    ceq(:, i) = x_m(1:4) - 1/2*(x_1(1:4) + x_2(1:4)) - h_k/8*(f_1-f_2);
    f_m = dynamics(t(i+1), x_m);
    ceq(:, i+1) = h_k/6 * (f_2 + 4*f_m + f_1) +x_1(1:4) - x_2(1:4);
end

ceq = reshape(ceq, [], 1);

end


