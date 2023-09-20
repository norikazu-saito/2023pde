function p1fem1(datafile, givenfunc)
    % read mesh data
    [deg, p, idp, t, idt, e, ide, d, idd] = getmesh(datafile);
    % coefficients 問題応じて定数を適切に設定すること
    coef_nu = 1.0; coef_c = 1.0; coef_kappa = 1.0;
    % coef_nu = 1.0; coef_c = 0.0; coef_kappa = 0.0;
    % assemble matrices and vectors
    [K, M, b] = matrix1(p, t, givenfunc);
    [R, rho] = robin1(p, e, coef_kappa, givenfunc);
    Aglobal = coef_nu*K + coef_c*M + coef_kappa*R;
    bglobal = b + rho;
    % solve the linear system
    [A, b, usol] = dirichlet1(p, d, Aglobal, bglobal, givenfunc);
    % 3D plot
    figure(5); plotP1fem(p, t, usol); saveas(5,"p1fem1.pdf");
    % output results
    F1 = fopen("p1fem1.res","w"); Z = [p;usol'];
    fprintf(F1,"%f %f %f\n",Z); fclose(F1);
end
