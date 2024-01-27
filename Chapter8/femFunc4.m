function z=femFunc4(id, x, y)
switch id
    case 0 % right-hand side funtion 
    z = zeros(size(x));
    case 1 % initial data
    %z = w;
    cx=1/4; cy=1/4; p=0.01; M=0.3; r=(x-cx).^2+(y-cy).^2;
    z=M*exp(-sqrt(r)/p);
    case 2 % b_1 
    z=sin(pi*x);
    case 3 % b_2 
    z=sin(pi*y);
    case 4 % exact solution 
    z = zeros(size(x)); 
end    
end