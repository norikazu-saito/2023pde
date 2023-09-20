function [K, M, C, kappa, bmax] = matrix1BT(p,t,givenfunc)
np = size(p,2); nt = size(t,2); kappa=100; bmax=-1;
K=sparse(np,np); M=sparse(np,np); C=sparse(np,np); force=zeros(np,1);
for l = 1:nt
  tlocal=t(1:3,l); x=p(1,tlocal); y=p(2,tlocal); [area,b,c]=P1grad(x,y);  
  % stiffness matrix
  Klocal=(b*b'+c*c')*area; K(tlocal,tlocal)=K(tlocal,tlocal)+ Klocal; 
  % lumping-mass matrix
  Mlocal = [4 0 0; 0 4 0; 0 0 4]/12*area; 
  M(tlocal,tlocal)=M(tlocal,tlocal)+Mlocal; 
  % convection matrix: Baba-Tabata upwinding
  xg=sum(x)/3; yg=sum(y)/3; 
  vecg=[xg;yg]; G=[vecg vecg vecg]; 
  bval=[givenfunc(2,xg,yg),givenfunc(3,xg,yg)];
  mx=(x+circshift(x,-1))/2; my=(y+circshift(y,-1))/2; 
  R=[mx;my]; S=[0 -1; 1 0]*(R-G); r=bval*S; rp=max(0,r); rm=max(0,-r);
  Clocal=[rp(1)+rm(3),-rm(1),-rp(3);-rp(1),rp(2)+rm(1),-rm(2);-rm(3),-rp(2),rp(3)+rm(2)];
  C(tlocal,tlocal)=C(tlocal,tlocal)+Clocal; 
  %perpendicular length
  xx=x-circshift(x,-1); yy=y-circshift(y,-1); 
  kappa=min(kappa,min(2*area*ones(size(xx))./sqrt(xx.^2+yy.^2)));
  %bmax
  bmax=max(bmax,norm(bval));
  %force term 
  %flocal=Mlocal*givenfunc(0,x,y)'; force(tlocal)=force(tlocal)+flocal; 
end
end
