function [hsize, l2, h1] = p1error(p, t, u, func)
    np = size(p,2); nt = size(t,2); sum1 = 0; sum2 = 0;
    for k = 1:nt
        % local coordinate
        tlocal = t(1:3,k); x = p(1,tlocal); y = p(2,tlocal);
        area = polyarea(x,y);
        % granularity parameter
        hh = sqrt((x-circshift(x,1)).^2 + (y-circshift(y,1)).^2);
        hsize = max(hh);
        % quadratue points & weights
        [w,Q] = gauss2D(x,y,area);
        % P1 function
        a = [ones(3,1),x',y'] \ u(tlocal); g = [ones(7,1),Q]*a;
        % function values
        f = (func(3,Q(:,1),Q(:,2)) - g).^2;
        dxf = (func(4,Q(:,1),Q(:,2)) - a(2)*ones(7,1)).^2;
        dyf = (func(5,Q(:,1),Q(:,2)) - a(3)*ones(7,1)).^2;
        % quadratue
        sum1 = sum1 + dot(f,w); sum2 = sum2 + dot(dxf+dyf,w);
    end
    l2 = sqrt(sum1); h1 = sqrt(sum2);  
end
