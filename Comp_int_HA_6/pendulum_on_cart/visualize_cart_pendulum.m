function visualize_cart_pendulum(x)
    l = length(x);
    x = reshape(x, 5, l/5);
    [~, n] = size(x);
    figure
    for i = 1:n
        clf

        xx0 = x(1,i);
        tta0 = x(2,i);

        hold on
        scatter(xx0,0,400,'k','s')
        plot([xx0, xx0 + sin(tta0)],[0, cos(tta0)],'k','LineWidth',2)
        scatter(xx0 + sin(tta0), cos(tta0),100,'k','filled')
        plot(x(1,1:i) + sin(x(2, 1:i)), cos(x(2,1:i)))
        xlim([-2.5 2.5])
        ylim([-2.5 2.5])

        pause(1e-1)
    end

end