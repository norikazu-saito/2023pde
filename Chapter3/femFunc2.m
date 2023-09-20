function z = femFunc2(id, x, y)
    w = sin(pi*x).*sin(pi*y);
    switch id
        case 0  % right-hand side function
            z = 2*pi^2*w;
        case 1  % Dirichlet g
            z = w;
        case 2  % Neumann r
            z = 0;
        case 3  % exact solution u
            z = w;
        case 4  % u_x
            z = pi*cos(pi*x).*sin(pi*y);
        case 5  % u_y 
            z = pi*sin(pi*x).*cos(pi*y);
    end
end
