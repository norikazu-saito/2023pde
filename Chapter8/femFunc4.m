function z=femFunc4(id, x, y)
w=cos(pi*x).*cos(pi*y)+1.2*ones(size(x)); 
switch id
    case 0 % right-hand side funtion 
    %z = 3*w + cos(x).*(cos(y).^2-sin(y).^2) + cos(y).*(cos(x).^2-sin(x).^2); 
    z = w; 
    case 1 % initial data
    %z = w;
    cx=1/4; cy=1/4; p=0.01; M=0.3; r=(x-cx).^2+(y-cy).^2;
    z=M*exp(-r/p);
    case 2 % b_1 
    %z = zeros(size(x)); 
    %z=(1-x).*x.*(1-2.*y);
    z=sin(pi*x);
    case 3 % b_2 
    %z = zeros(size(x));
    %z=-(1-y).*y.*(1-2.*x);
    z=sin(pi*y);
    case 4 % exact solution 
    z = w; 
end    
end