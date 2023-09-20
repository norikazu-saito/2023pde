function [K, M, force] = matrix1(p, t, givenfunc)
    np = size(p,2); nt = size(t,2);
    K = sparse(np,np); M = sparse(np,np); force = zeros(np,1);
    for l = 1:nt
        tlocal = t(1:3,l); x = p(1,tlocal); y = p(2,tlocal);
        [area,b,c] = P1grad(x,y);
        % stiffness matrix
        Klocal = (b*b'+c*c')*area;
        K(tlocal,tlocal) = K(tlocal,tlocal) + Klocal;
        % mass matrix
        Mlocal = [2 1 1; 1 2 1; 1 1 2]/12*area;
        M(tlocal,tlocal) = M(tlocal,tlocal) + Mlocal;
        % force term
        flocal = Mlocal*givenfunc(0,x,y)';
        force(tlocal) = force(tlocal) + flocal;
    end
end
