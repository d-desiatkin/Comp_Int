function [fun_val] = MPSO(fun, x_dim, n_particles, inert_w, part_w, global_w, blo, bup)
%MPSO Summary of this function goes here
%   Detailed explanation goes here

x = blo + (bup-blo)*rand(x_dim, n_particles);
p = x;
g = p(:,1);
for s = 2:n_particles
    if(fun(p(:,s)) < fun(g))
        g = p(:,s);
    end
end

v = -abs(bup-blo) + abs(bup-blo)*rand(x_dim, n_particles);

n = 1;

while 1
    
    
    for s = 1:n_particles
        rp = rand(x_dim,1);
        rg = rand(x_dim,1);
        v(:,s) = inert_w * v(:,s) + part_w * rp .* (p(:,s)-x(:,s)) + ...
            global_w*rg.*(g-x(:,s));
        x(:,s) = x(:,s) + v(:,s);
        if fun(x(:,s)) < fun(p(:,s))
            p(:,s) = x(:,s);
            if fun(p(:,s)) < fun(g)
                g = p(:,s);
            end
        end 
    end
    
    if n == 1000
        break;
    end
    
    n = n + 1;
    
end

fun_val = g;

end

