function dxdt = abg_dnmcs(t,x)
    dxdt = zeros(2,1);
    
    s = x(1);
    s_d = x(2);
    
    dxdt(1) = s_d;
    dxdt(2) = -( gamma_fcn(s) + beta_fcn(s)*s_d^2 )/alfa_fcn(s);
end