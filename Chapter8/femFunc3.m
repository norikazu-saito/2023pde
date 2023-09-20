function z=femFunc3(id, x, y) 
switch id
    case 0 % right-hand side function
    z = ones(size(x)); 
    case 1 % Dirichlet g 
    z = zeros(size(x));
    case 2 % Neumann r 
    z = zeros(size(x)); 
    case 3 % exact solution u 
    z = zeros(size(x));
    case 4 % u_x 
    z = zeros(size(x));
    case 5 % u_y 
    z = zeros(size(x));
end    
end