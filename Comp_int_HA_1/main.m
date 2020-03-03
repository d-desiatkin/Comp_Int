pkg load geometry
links = [1.075, 1.280]
joints = [pi/9, -pi/9, pi/9, pi/9, -pi/9, pi/9]
init_joints = [1, 1, 1, 1, 1, 1]
con_a = [-pi, -pi/3, -5/6*pi, -2*pi, -2/3*pi, -2*pi];

con_b = [pi, pi/3, 5/6*pi, 2*pi, 2/3*pi/9, 2*pi];

H_base = [1 0 0 0;
          0 1 0 0;
          0 0 1 0.670;
          0 0 0 1];
       
H_target = H_base * Rz(joints(1)) * Tx(0.312) * Ry(joints(2)) * Tz(links(1)) * ...
    Ry(joints(3)) * Tz(0.225) * Tx(links(2)) * Rx(joints(4)) * ...
    Ry(joints(5)) * Rx(joints(6)) * Tx(0.215)
    
     
options = optimset('Display', 'iter', 'GradObj', 'off', 'MaxIter', 400);
[j,cost] = fmincon(@(q)(CostFun(q, H_base, links, H_target)), init_joints, options)

H_pred = H_base * Rz(j(1)) * Tx(0.312) * Ry(j(2)) * ...
         Tz(links(1)) * Ry(j(3)) * Tz(0.225) * Tx(links(2)) * ...
         Rx(j(4)) * Ry(j(5)) * Rx(j(6)) * ...
         Tx(0.215)