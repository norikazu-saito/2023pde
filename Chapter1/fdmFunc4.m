function z=fdmFunc4(id, x, t)
switch id
    case 0 % initial value 
        z= min(x,ones(size(x))-x);
    case 1 % source term 
        z = zeros(size(x)); 
    case 2 % exact solution
        z = zeros(size(x)); 
end    
end