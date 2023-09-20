function [K, M, B, force] = matrixStab1(p,t,convect,coef_c,delta,givenfunc)
np = size(p,2); nt = size(t,2); K = sparse(np,np); M = sparse(np,np); 
B = sparse(np,np); force = zeros(np,1); bmax=norm(convect);
%    
for l = 1:nt
  tlocal=t(1:3,l);x=p(1,tlocal);y=p(2,tlocal);[area,b,c]=P1grad(x,y);  
  % granularity parameter
  hh=sqrt((x-circshift(x,1)).^2 + (y-circshift(y,1)).^2); hsize=max(hh);
  % stiffness matrix
  Klocal=(b*b'+c*c')*area; K(tlocal,tlocal)=K(tlocal,tlocal)+ Klocal; 
  % mass matrix
  Mlocal = [2 1 1; 1 2 1; 1 1 2]/12*area; 
  M(tlocal,tlocal)=M(tlocal,tlocal)+Mlocal; 
  % convection matrix 
  Bxlocal=[b b b]'*area/3; Bylocal=[c c c]'*area/3;
  B(tlocal,tlocal)=B(tlocal,tlocal) + convect(1)*Bxlocal + convect(2)*Bylocal;
  %stabilization
  tau=delta*hsize/bmax;
  Bstab=convect(1)*coef_c*Bxlocal+convect(2)*coef_c*Bylocal;
  Bstab=Bstab+convect(1)^2*b*b'*area+convect(1)*convect(2)*b*c'*area+convect(1)*convect(2)*c*b'*area+convect(2)^2*c*c'*area;
  B(tlocal,tlocal)=B(tlocal,tlocal) + tau*Bstab;
  % force term 
  flocal=Mlocal*givenfunc(0,x,y)'; force(tlocal)=force(tlocal)+flocal;
  flocal=(convect(1)*Bxlocal+convect(2)*Bylocal)*givenfunc(0,x,y)'; 
  force(tlocal)=force(tlocal)+tau*flocal;
end
end
