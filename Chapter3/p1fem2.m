function [hsize, l2, h1, linf] = p1fem2(datafile, givenfunc, id)
    %% if id == 0, then return and plot usol
    %% if id > 0, then, in addition above, return errors using p1error.m
    % get mesh data
    [deg, p, idp, t, idt, e, ide, d, idd] = getmesh(datafile);
    % coefficients
    %% coef_nu = 1.0; coef_c = 1.0; coef_kappa = 1.0;
    coef_nu = 1.0; coef_c = 0.0; coef_kappa = 0.0;
    % assemble matrices and vectors
    [K, M, b] = matrix1(p, t, givenfunc);
    [R, rho] = robin1(p, e, coef_kappa, givenfunc);
    Aglobal = coef_nu*K + coef_c*M + coef_kappa*R;
    bglobal = b + rho;
    % solve the linear system
    [A, b, usol] = dirichlet1(p, d, Aglobal, bglobal, givenfunc);
    % 3D plot
    figure(5); plotP1fem(p,t,usol); saveas(5,"p1fem2.pdf");
    % output results
    F1 = fopen("p1fem2.res","w"); Z = [p;usol']; fprintf(F1,"%f %f %f\n",Z);
    fclose(F1);
    % error observation
    if id == 0
        hsize = 1; l2 = 1; h1 = 1; linf = 1;
    else
        [hsize, l2, h1] = p1error(p, t, usol, givenfunc);
        linf = norm(usol'-givenfunc(3,p(1,:),p(2,:)),inf);
    end
end
