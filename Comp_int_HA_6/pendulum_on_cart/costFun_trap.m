function [J] = costFun_trap(t, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
l = length(x);
x = reshape(x, 5, l/5);
J = 0;
for i=1:l/5-1
    h_k = (t(i+1)-t(i));
    
    J = J + h_k/2*(x(5, i)^2 + x(5, i+1)^2);

end
end

