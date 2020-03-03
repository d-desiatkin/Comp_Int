function [c,ceq] = mycon(t, x)
%MYCON Summary of this function goes here
%   Detailed explanation goes here
c = [];  % Compute nonlinear inequalities at x.
ceq = [];   % Compute nonlinear equalities at x.

l = length(x);
x = reshape(x, 5, l/5);
[~, n] =  size(x);

ceq = zeros(4, n-1);

for i=1:n-1
    x_1 = x(:, i);
    x_2 = x(:, i+1);
    f_1 = dynamics(t(i), x_1);
    f_2 = dynamics(t(i+1), x_2);
    ceq(:,i) = (t(i+1)-t(i))/2 * (f_2 + f_1) +x_1(1:4) - x_2(1:4);
end

ceq = reshape(ceq, [], 1);

end


