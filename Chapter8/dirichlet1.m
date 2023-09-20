function [A, b, u] = dirichlet1(p, d, A, b, givenfunc)
    np = size(p,2); u = zeros(np,1);
    % boundary and interior nodes
    fixed = unique(d); free = setdiff([1:np],fixed);
    % set boundary value
    g = zeros(np,1); x = p(1,fixed); y = p(2,fixed);
    g(fixed) = givenfunc(1,x,y);
    % re-allocate A and b
    b = b(free) - A(free,fixed)*g(fixed); A = A(free,free);
    % solving Au = b
    u(fixed) = g(fixed); u(free) = A\b; 
end
