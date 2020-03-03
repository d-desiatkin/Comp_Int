function dxdt = dynamics(t,x)
    
    dxdt = zeros(4,1);
    
    Ks = 50000;
    Ku = 200000;
    Ms = 320;
    Mu = 40;
    Bs = 2000;
    
    M = [0 1 0 -1
        -Ks/Ms -Bs/Ms 0 Bs/Ms
         0 0 0 1
         Ks/Mu Bs/Mu -Ku/Mu -Bs/Mu];
    
    dxdt = M*x;
end