function [R, rho] = robin1(p, e, coef_kappa, givenfunc)
    np = size(p,2); ne = size(e,2);
    R = sparse(np,np); rho = zeros(np,1);
    for k = 1:ne
        elocal = e(1:2,k); x = p(1,elocal); y = p(2,elocal);
        length = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        % for mass matrix
        Rlocal = [2 1; 1 2]*length/6;
        R(elocal,elocal) = R(elocal,elocal) + Rlocal;
        % for right-hand side vector
        tmp = coef_kappa*givenfunc(1,x,y) + givenfunc(2,x,y);
        rlocal = Rlocal*tmp'; rho(elocal) = rho(elocal) + rlocal;
    end
end
