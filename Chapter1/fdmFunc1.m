function z = fdmFunc1(id, x, t)
    switch id
        case 0  % initial value
            z = x.*(sin(3*pi*x)).^2;
        case 1  % source term
            z = zeros(size(x));
        case 2  % exact solution
            z = zeros(size(x));
    end
end
