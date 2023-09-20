function z=fdmFunc3(id, x, t)
switch id
    case 0 % initial value 
        I=find(x<=0.3); J=find(x>0.3 & x <=0.6); K=find(x>0.6); 
        z=zeros(size(x)); z(I)=0.3; z(J)=-2*(x(J)-0.6); z(K)=0.6;
    case 1 % source term 
        z = zeros(size(x)); 
    case 2 % exact solution
        z = zeros(size(x)); 

end    
end