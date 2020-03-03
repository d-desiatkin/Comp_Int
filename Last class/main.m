

A = [-1, 0;
     0, -2];


cvx_begin sdp

    variable P(2,2) symmetric
    subject to
        %P > 0;
        A'*P + P*A < -eye(2);

cvx_end


P