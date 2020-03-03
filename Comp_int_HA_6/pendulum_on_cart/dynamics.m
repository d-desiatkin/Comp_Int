function dxdt = dynamics(t,x)
    dxdt = zeros(4,1);
    
    M = M_mtrx(x(1:2));
    n = n_vctr(x(1:2), x(3:4));
    
    dxdt(1:2) = x(3:4);
    dxdt(3:4) = M\([1, 0]'*x(5) - n);
end

