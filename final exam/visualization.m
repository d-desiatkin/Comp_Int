function visualization(x,l)
    figure()
    for i = 1:length(x)
        clf
        q1=x(i,1);
        q2=x(i,2);
        
        DH=[q1,0,l(1),0;
            q2,0,l(2),0;];
        
        H=eye(4);
        x1=[0];
        y1=[0];
        z1=[0];
        n=size(DH);

for i=1:n(1)
    H=H*Rz(DH(i,1))*Tz(DH(i,2))*Tx(DH(i,3))*Rx(DH(i,4));
    x1(i+1)=H(1,4);
    y1(i+1)=H(2,4);
    z1(i+1)=H(3,4);
end

        hold on
        plot(x1,y1,'k')
        xlabel('x')
        ylabel('y')
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        
        pause(1e-4)
    end
end

