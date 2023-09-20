function z=fdmFunc5(id, x, t)
switch id
    case 0 % initial value 
        z = x .^3 .* (1-x);
    case 1 % source term 
        z = exp(t) .* (-x.^4 + x.^3 + 12*x.^2 -6*x);
    case 2 % exact solution
        z = exp(t) .* x .^3 .* (1-x); 
end    
end