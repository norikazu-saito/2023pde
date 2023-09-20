function heatRate(N, Tmax, theta, lam, imax, givenfunc)
    % preliminaries
    hv = []; e1 = []; hv1 = []; e2 = []; N1 = N;
    % convergence
    for i = 1:imax
        [h, err1, err2] = heat(N1, theta, lam, Tmax, givenfunc, -1);
        hv = [hv;h]; e1 = [e1;err1]; e2 = [e2;err2];
        if i == 1
            sl1 = 0.1*err2; hv1 = [hv1;sl1];
        else
            sl1 = sl1/4; hv1 = [hv1;sl1];
        end
        % granulating
        N1 = N1*2;
    end
    % log plot
    figure(1); loglog(hv,e1,'-xb','LineWidth',2); hold on;
    loglog(hv,e2,':sr','LineWidth',2); loglog(hv,hv1,'-k','LineWidth',1);
    hold off;
    % decoration
    xlabel('h'); ylabel('ERR'); grid on; pbaspect([1 2 1]);
    legend('L^{inf} err','L^2 err','slope 2','Location','SouthEast');
    saveas(1,'heatRate.pdf');
end
