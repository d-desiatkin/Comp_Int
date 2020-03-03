function [dxdt, ftraj, fu] = dynamics(t, x, K, x_0)
    
    if (nargout == 3) | (nargout == 2)
        [n, d] = size(x);
        fu = zeros(n,round(d/2));
        ftraj = zeros(n,round(d/2));
        dxdt = 0;
        for i=1:n
            traj_p = [sin(t(i,1)) * pi/180*15 + x_0(1); cos(t(i,1))*pi/180*15 + x_0(2);];
            u_p = -K*(traj_p - x(i,1:2)');
            fu(i,:) = u_p';
            ftraj(i,:) = traj_p';
        end
        
    end
    
    if(nargout == 1)
        traj_p = [sin(t) * pi/180*20 + x_0(1); cos(t)*pi/180*20 + x_0(2);];
        dxdt = zeros(4,1);
        u_p = -K*(x(1:2) - traj_p);
        M = M_mtrx(x(1:2));
        n = n_vctr(x(1:2), x(3:4));
        dxdt(1:2) = x(3:4);
        dxdt(3:4) = M\(u_p-n);
    end
    
end