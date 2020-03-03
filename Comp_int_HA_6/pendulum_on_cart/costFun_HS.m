function [J] = costFun_HS(t, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
l = length(x);
x = reshape(x, 5, l/5);
J = 0;
for i=1:2:l/5-2
    h_k = (t(i+2)-t(i));
    
    J = J + h_k/6*(x(5, i)^2 + 4*x(5, i+1)^2 + x(5, i+2)^2);

end
end

