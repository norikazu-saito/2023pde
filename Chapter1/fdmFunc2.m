function z = fdmFunc2(id, x, t)
    switch id
        case 0  % initial value
            z = sin(pi*x);
        case 1  % source term
            z = zeros(size(x));
        case 2  % exact solution
            z = exp(-t*pi^2)*sin(pi*x);
    end
end
