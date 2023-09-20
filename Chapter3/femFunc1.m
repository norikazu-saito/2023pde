function z = femFunc1(id, x, y)
    a = 0.5; b = 0.5; lam = 5; r = (x-a).^2 + (y-b).^2; w = exp(-lam*r);
    switch id
        case 0  % right-hand side function
            z = w .* (1 + 4*lam - 4*lam^2.*r);
        case 1  % Dirichlet g
            z = w;
        case 2  % Neumann r
            z = 2 * lam * (y-b) .* w;
        case 3  % exact solution u
            z = w;
        case 4  % exact solution u_x
            z = -2 * lam * (x-a) .* w;
        case 5  % exact solution u_y
            z = -2 * lam * (y-b) .* w;
    end
end
