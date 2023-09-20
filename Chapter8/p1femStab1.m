function [hsize,l2,h1,linf] = p1femStab1(datafile,givenfunc,coef_nu,coef_c,delta,id)
%% if id==0, then return and plot usol 
%% if id>0,  then, in addition above, retuen errors using p1error.m 
% get mesh data
[deg,p,idp,t,idt,e,ide,d,idd] = getmesh(datafile);
% convection term
convect=[5.0,5.0]';
% assemble matrices and vectors
[K, M, B, b] = matrixStab1(p,t,convect,coef_c,delta,givenfunc);
Aglobal = coef_nu*K + B + coef_c*M;
% solve the linear system
[A,b,usol] = dirichlet1(p,d,Aglobal,b,givenfunc); 
% 3D plot
figure(5); plotP1fem(p,t,usol); saveas(5,"p1femStab1.pdf");
% output resuts
F1=fopen("p1femStab1.res","w"); Z=[p;usol']; 
fprintf(F1,"%f %f %f\n",Z); 
fclose(F1);
% error observtion 
if id==0 
    hsize=1; l2=1; h1=1; linf=1;  
else 
    [hsize,l2,h1] = p1error(p,t,usol,givenfunc);
    linf=norm(usol'-givenfunc(3,p(1,:),p(2,:)),inf);
end
end